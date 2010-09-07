/* -*- Mode: C -*- */

#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"

#include "ppport.h"

#undef die

#include <gaul.h>
#include "mine.h"

static HV *population_stash;

/* user data AV slots */

struct my_data {
    SV *wrapper;
    int chromosome_bytes;
    SV **chromosome_svs;
    byte *bytes;
    unsigned int max_bytes;
    SV *evaluate_cb;
    SV *adapt_cb;
    SV *stop_cb;
    SV *select_one_cb;
    SV *select_two_cb;
    SV *replace_cb;
    SV *seed_cb;
    SV *crossover_cb;
    SV *mutate_cb;
};

#define pop_data ((struct my_data *)((pop)->data))

static SV*
get_chromosome_sv(population *pop, int i) {
    SV **svs = pop_data->chromosome_svs;
    SV *sv = svs[i];
    if (SvREFCNT(sv) > 1) {
        SV *alt = svs[pop->num_chromosomes];
        if (SvREFCNT(alt) > 1) {
            fprintf(stderr, "(sv realloc)"); fflush(stderr);
            sv_2mortal(sv);
            return pop_data->chromosome_svs[i] = newSV(pop_data->chromosome_bytes);
        }
        else {
            svs[pop->num_chromosomes] = sv;
            svs[i] = alt;
            return alt;
        }
    }
    return sv;
}

enum slots_ {
    av_slot_ptr = 0,
    av_slot_evaluate,
    av_slot_seed,
    av_slot_adapt,
    av_slot_select_one,
    av_slot_select_two,
    av_slot_mutate,
    av_slot_crossover,
    av_slot_replace,
    av_slot_size
};

struct enum_table {
    char *name;
    int value;
};

typedef void (*cb_C)(void);

struct cb_table {
    char *name;
    cb_C cb;
};

static boolean
my_ga_chromosome_bitstring_allocate(population *pop, entity *embryo)  {
    int i;		/* Loop variable over all chromosomes */
    int bytes;

    if (!pop) die("Null pointer to population structure passed.");
    if (!embryo) die("Null pointer to entity structure passed.");

    if (embryo->chromosome != NULL)
        die("This entity already contains chromosomes.");

    embryo->chromosome = s_malloc(pop->num_chromosomes * sizeof(byte *));
    
    bytes = (pop->len_chromosomes + 7) / 8;
    for (i = 0; i < pop->num_chromosomes; i++) {
        byte *chromo = (byte *)s_malloc(bytes);
        if (!chromo) die("Unable to allocate bitstring");
        if (bytes) chromo[bytes - 1] = 0; /* clear extra bits in the last byte */
        embryo->chromosome[i] = chromo;
    }
    return TRUE;
}

static unsigned int
my_ga_chromosome_bitstring_to_bytes(const population *pop,
                                    entity *joe,
                                    byte **bytes, unsigned int *max_bytes) {
    int num_bytes;	/* Actual size of genes. */
    int i;		/* Loop variable over all chromosomes */

    if (!pop) die("Null pointer to population structure passed.");
    if (!joe) die("Null pointer to entity structure passed.");


    num_bytes = ga_bit_sizeof(pop->len_chromosomes * pop->num_chromosomes);

    if (num_bytes > *max_bytes) {
        *max_bytes = num_bytes;
        *bytes = s_realloc(*bytes, *max_bytes);
    }

    (*bytes)[(num_bytes > 0) ? num_bytes - 1 : 0] = 0;

    if (!joe->chromosome) {
        return 0;
    }

    for(i=0; i<pop->num_chromosomes; i++) {
        ga_bit_copy( *bytes, joe->chromosome[i],
                     i*pop->len_chromosomes, 0,
                     pop->len_chromosomes );
    }

    return num_bytes;
}

static SV *
entity_to_sv(population *pop, entity *joe) {
    SV *sv;
    /* fprintf(stderr, "bytes: %p, max_bytes: %d\n", pop_data->bytes, pop_data->max_bytes); fflush(stderr); */
    int len = pop->chromosome_to_bytes(pop, joe, &(pop_data->bytes), &(pop_data->max_bytes));
    return newSVpvn(pop_data->bytes, len);
}

static int
sv_to_enum(SV *sv, struct enum_table *table, const char *enum_name, int def) {
    if (SvOK(sv)) {
        STRLEN len;
        const char *pv = SvPV_const(sv, len);
        for (; table->name; table++) {
            if (!strcmp(table->name, pv))
                return table->value;
        }
        croak("bad value for %s enumeration", enum_name);
    }
    return def;
}

static cb_C
sv_to_cb(SV *sv, struct cb_table *table, const char *cb_name, SV **slot, const char *def) {
    const char *name;
    SV *perl_cb = NULL;
    if (SvOK(sv)) {
        if (SvROK(sv) && SvTYPE(SvRV(sv)) == SVt_PVCV) {
            name = "_perl";
            perl_cb = sv;
        }
        else {
            name = SvPV_nolen_const(sv);
            if (name[0] == '_') croak("bad value for %s callback", cb_name);
        }
    }
    else name = def;

    if (name) {
        for (; table->name; table++) {
            if (!strcmp(table->name, name)) {
                if (perl_cb) *slot = SvREFCNT_inc(perl_cb);
                return table->cb;
            }
        }
        croak("bad value for %s callback", cb_name);
    }
    return NULL;
}

#include "enum-tables.h"

#include "cb.h"
#include "cb-tables.h"

static SV *
delete_opt_(HV *opts, const char *name, STRLEN name_len, int required) {
    SV *v = hv_delete(opts, name, name_len, 0);
    if (!v) {
        if (required)
            croak("required argument %s missing", name);
        //_exit(4);
        return &PL_sv_undef;
    }
    if (required && !SvOK(v))
        croak("required argument %s in undef", name);
    return v;
}

#define delete_opt_required(name) delete_opt_(opts, name, strlen(name), 1)
#define delete_opt(name) delete_opt_(opts, name, strlen(name), 0)

static int
delete_opt_int_(HV *opts, const char *name, STRLEN name_len, int def) {
    SV *sv = delete_opt(name);
    return (SvOK(sv) ? SvIV(sv) : def);
}

#define delete_opt_int(name, def) delete_opt_int_(opts, name, strlen(name), def)
#define delete_opt_int_required(name) SvIVx(delete_opt_required(name))

typedef GAgeneration_hook GAstop;

#define delete_opt_cb(name, def) (GA ## name)sv_to_cb(delete_opt(#name), table_for_ga_ ## name, #name, &(pop_data->name ## _cb), def)

#define delete_opt_cb_with_key(table, key, def) (GA ## key)sv_to_cb(delete_opt(#key), table_for_ga_ ## table, #key, &(pop_data->key ## _cb), def)

MODULE = AI::GAUL		PACKAGE = AI::GAUL		

population *
genesis_any(klass, opts)
    SV *klass
    HV *opts
PREINIT:
    population *pop;
    int num_chromo;
    int len_chromo;
    int population_size;
    int type_chromo;
    SV *sv;
    int i;
    struct my_data *data;
CODE:
{
    population_size = delete_opt_int_required("population_size");
    len_chromo = delete_opt_int_required("len_chromo");
    num_chromo = delete_opt_int("num_chromo", 1);

    pop = ga_population_new(population_size, num_chromo, len_chromo);
    if (!pop)
        croak("Unable to allocate population object");

    Newxz(data, 1, struct my_data);
    pop->data = data;

    pop->iteration_hook = NULL;
    pop->data_destructor = NULL;
    pop->data_ref_incrementor = NULL;

    pop->evaluate = delete_opt_cb(evaluate, NULL);
    pop->adapt = delete_opt_cb(adapt, NULL);
    pop->select_one = delete_opt_cb(select_one, "sus");
    pop->select_two = delete_opt_cb(select_two, "sus");
    pop->replace = delete_opt_cb(replace, NULL);

    pop->generation_hook = delete_opt_cb(stop, NULL);

    sv = delete_opt("scheme");
    if (SvOK(sv))
        ga_population_set_scheme(pop, sv_to_enum(sv, table_for_ga_scheme_type_t, "scheme", 0));

    sv = delete_opt("elitism");
    if (SvOK(sv))
        ga_population_set_elitism(pop, sv_to_enum(sv, table_for_ga_elitism_type_t, "elitism", 0));

    sv = delete_opt("crossover_prob");
    if (SvOK(sv))
        ga_population_set_crossover(pop, SvNV(sv));

    sv = delete_opt("mutation_prob");
    if (SvOK(sv))
        ga_population_set_mutation(pop, SvNV(sv));

    sv = delete_opt("migration_prob");
    if (SvOK(sv))
        ga_population_set_mutation(pop, SvNV(sv));

    sv = delete_opt("allele_mutation_prob");
    if (SvOK(sv))
        ga_population_set_allele_mutation_prob(pop, SvNV(sv));

    sv = delete_opt("type_chromo");
    switch (sv_to_enum(sv, table_for_chromo_type_t, "type_chromo", GA_CHROMO_TYPE_BITSTRING)) {
    case GA_CHROMO_TYPE_BITSTRING:
        pop->chromosome_constructor = my_ga_chromosome_bitstring_allocate;
        pop->chromosome_destructor = ga_chromosome_bitstring_deallocate;
        pop->chromosome_replicate = ga_chromosome_bitstring_replicate;
        pop->chromosome_to_bytes = my_ga_chromosome_bitstring_to_bytes;
        pop->chromosome_from_bytes = ga_chromosome_bitstring_from_bytes;
        pop->chromosome_to_string = ga_chromosome_bitstring_to_string;
        pop->seed = delete_opt_cb_with_key(seed_bitstring, seed, "random");
        pop->crossover = delete_opt_cb_with_key(crossover_bitstring, crossover, "mixing");
        pop->mutate = delete_opt_cb_with_key(mutate_bitstring, mutate, "multipoint");
        data->chromosome_bytes = (len_chromo + 7) / 8;
        break;
    case GA_CHROMO_TYPE_BOOLEAN:
        pop->chromosome_constructor = ga_chromosome_boolean_allocate;
        pop->chromosome_destructor = ga_chromosome_boolean_deallocate;
        pop->chromosome_replicate = ga_chromosome_boolean_replicate;
        pop->chromosome_to_bytes = ga_chromosome_boolean_to_bytes;
        pop->chromosome_from_bytes = ga_chromosome_boolean_from_bytes;
        pop->chromosome_to_string = ga_chromosome_boolean_to_string;
        pop->seed = delete_opt_cb_with_key(seed_boolean, seed, "random");
        pop->crossover = delete_opt_cb_with_key(crossover_boolean, crossover, "mixing");
        pop->mutate = delete_opt_cb_with_key(mutate_boolean, mutate, "multipoint");
        data->chromosome_bytes = len_chromo * sizeof(boolean);
        break;
    case GA_CHROMO_TYPE_CHAR:
        pop->chromosome_constructor = ga_chromosome_char_allocate;
        pop->chromosome_destructor = ga_chromosome_char_deallocate;
        pop->chromosome_replicate = ga_chromosome_char_replicate;
        pop->chromosome_to_bytes = ga_chromosome_char_to_bytes;
        pop->chromosome_from_bytes = ga_chromosome_char_from_bytes;
        pop->chromosome_to_string = ga_chromosome_char_to_string;
        pop->seed = delete_opt_cb_with_key(seed_char, seed, "random");
        pop->crossover = delete_opt_cb_with_key(crossover_char, crossover, "mixing");
        pop->mutate = delete_opt_cb_with_key(mutate_char, mutate, "multipoint");
        data->chromosome_bytes = len_chromo * sizeof(char);
        break;
    case GA_CHROMO_TYPE_DOUBLE:
        pop->chromosome_constructor = ga_chromosome_double_allocate;
        pop->chromosome_destructor = ga_chromosome_double_deallocate;
        pop->chromosome_replicate = ga_chromosome_double_replicate;
        pop->chromosome_to_bytes = ga_chromosome_double_to_bytes;
        pop->chromosome_from_bytes = ga_chromosome_double_from_bytes;
        pop->chromosome_to_string = ga_chromosome_double_to_string;
        pop->seed = delete_opt_cb_with_key(seed_double, seed, "random");
        pop->crossover = delete_opt_cb_with_key(crossover_double, crossover, "mixing");
        pop->mutate = delete_opt_cb_with_key(mutate_double, mutate, "multipoint");
        data->chromosome_bytes = len_chromo * sizeof(double);
        sv = delete_opt("allele_min");
        if (SvOK(sv))
            ga_population_set_allele_min_double(pop, SvNV(sv));
        sv = delete_opt("allele_max");
        if (SvOK(sv))
            ga_population_set_allele_max_double(pop, SvNV(sv));

        break;
    case GA_CHROMO_TYPE_INTEGER:
        pop->chromosome_constructor = ga_chromosome_integer_allocate;
        pop->chromosome_destructor = ga_chromosome_integer_deallocate;
        pop->chromosome_replicate = ga_chromosome_integer_replicate;
        pop->chromosome_to_bytes = ga_chromosome_integer_to_bytes;
        pop->chromosome_from_bytes = ga_chromosome_integer_from_bytes;
        pop->chromosome_to_string = ga_chromosome_integer_to_string;
        pop->seed = delete_opt_cb_with_key(seed_integer, seed, "random");
        pop->crossover = delete_opt_cb_with_key(crossover_integer, crossover, "mixing");
        pop->mutate = delete_opt_cb_with_key(mutate_integer, mutate, "multipoint");
        data->chromosome_bytes = len_chromo * sizeof(int);
        sv = delete_opt("allele_min");
        if (SvOK(sv))
            ga_population_set_allele_min_integer(pop, SvIV(sv));

        sv = delete_opt("allele_max");
        if (SvOK(sv))
            ga_population_set_allele_max_integer(pop, SvIV(sv));
        
        break;
    default:
        croak("bad chromo type");
    }

    Newxz(data->chromosome_svs, num_chromo + 1, SV*);
    for (i = 0; i <= num_chromo; i++)
        data->chromosome_svs[i] = newSV(data->chromosome_bytes);

    ga_population_seed(pop);

    RETVAL = pop;
}
OUTPUT:
    RETVAL


MODULE = AI::GAUL		PACKAGE = AI::GAUL    PREFIX = ga_

int
ga_evolution(pop, max_generations)
    population *pop
    int max_generations

int
ga_evolution_steady_state(pop, max_generations)
    population *pop
    int max_generations

void
ga_dump(pop)
    population *pop
CODE:
    ga_population_dump(pop);

SV *
ga_fitness_mean(pop)
    population *pop
PREINIT:
    double average;
CODE:
    if (ga_fitness_mean(pop, &average))
        RETVAL = newSVnv(average);
    else
        RETVAL = &PL_sv_undef;
OUTPUT:
    RETVAL

void
ga_fitness_stats(pop)
    population *pop
PREINIT:        
    double max, min, mean, median, variance, stddev, kurtosis, skew;
PPCODE:
    if (ga_fitness_stats(pop, &max, &min, &mean, &median, &variance, &stddev, &kurtosis, &skew)) {
        mXPUSHn(max);
        mXPUSHn(min);
        mXPUSHn(mean);
        mXPUSHn(median);
        mXPUSHn(variance);
        mXPUSHn(stddev);
        mXPUSHn(kurtosis);
        mXPUSHn(skew);
        XSRETURN(8);
    }
    else
        XSRETURN(0);

int
ga_id_from_rank(pop, rank)
    population *pop
    int rank
CODE:
    RETVAL = ga_get_entity_id_from_rank(pop, rank);
OUTPUT:
    RETVAL

void
ga_chromosomes_from_id(pop, id)
    population *pop
    int id
PREINIT:
    entity *joe;
PPCODE:
    if (joe = ga_get_entity_from_id(pop, id)) {
        int i;
        byte **chromo = (byte**)joe->chromosome;
        int bytes = pop_data->chromosome_bytes;
        for (i = 0; i < pop->num_chromosomes; i++)
            XPUSHs(chromo ? sv_2mortal(newSVpvn(chromo[i], bytes)) : &PL_sv_undef);
        XSRETURN(i);
    }
    else XSRETURN(0);

void
ga_chromosomes_from_rank(pop, rank)
    population *pop
    int rank
PREINIT:
    entity *joe;
PPCODE:
    if (joe = ga_get_entity_from_rank(pop, rank)) {
        int i;
        byte **chromo = (byte**)joe->chromosome;
        int bytes = pop_data->chromosome_bytes;
        for (i = 0; i < pop->num_chromosomes; i++)
            XPUSHs(chromo ? sv_2mortal(newSVpvn(chromo[i], bytes)) : &PL_sv_undef);
        XSRETURN(i);
    }
    else XSRETURN(0);

SV *
ga_fitness_from_rank(pop, rank)
    population *pop
    int rank
PREINIT:
    entity *joe;
CODE:
    if (joe = ga_get_entity_from_rank(pop, rank)) {
        RETVAL = newSVnv(joe->fitness);
    }
    else RETVAL = &PL_sv_undef;
OUTPUT:
    RETVAL

SV *
ga_fitness_from_id(pop, id)
    population *pop
    int id
PREINIT:
    entity *joe;
CODE:
    if (joe = ga_get_entity_from_id(pop, id)) {
        RETVAL = newSVnv(joe->fitness);
    }
    else RETVAL = &PL_sv_undef;
OUTPUT:
    RETVAL

int
ga_len_chromo(pop)
    population *pop
CODE:
    RETVAL = pop->len_chromosomes;
OUTPUT:
    RETVAL

int
ga_num_chromo(pop)
    population *pop
CODE:
    RETVAL = pop->num_chromosomes;
OUTPUT:
    RETVAL


void
DESTROY(pop)
    population *pop
PREINIT:
    struct my_data *data;
    SV **svs;
CODE:
    if (pop) {
        if (data = pop_data) {
            if (svs = data->chromosome_svs) {
                int i;
                for (i= 0; i <= pop->num_chromosomes; i++) sv_2mortal(svs[i]);
                Safefree(svs);
            }
            if (data->max_bytes) s_free(data->bytes);
            if (data->evaluate_cb) SvREFCNT_dec(data->evaluate_cb);
            if (data->adapt_cb) SvREFCNT_dec(data->adapt_cb);
            if (data->stop_cb) SvREFCNT_dec(data->stop_cb);
            if (data->select_one_cb) SvREFCNT_dec(data->select_one_cb);
            if (data->select_two_cb) SvREFCNT_dec(data->select_two_cb);
            if (data->replace_cb) SvREFCNT_dec(data->replace_cb);
            if (data->seed_cb) SvREFCNT_dec(data->seed_cb);
            if (data->mutate_cb) SvREFCNT_dec(data->mutate_cb);
            if (data->crossover_cb) SvREFCNT_dec(data->replace_cb);
            pop->data = NULL;
        }
        ga_extinction(pop);
        sv_setiv(SvRV(ST(0)), 0);
    }

void
ga_random_seed(seed)
    int seed


BOOT:
    ga_init_openmp();
    random_init();
    population_stash = gv_stashpv("AI::GAUL::Population", 1);


