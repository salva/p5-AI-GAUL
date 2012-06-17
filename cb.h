
static boolean
ga_evaluate__perl(population *pop, entity *joe) {
    dSP;
    SV *ret_sv;
    boolean ret;
    int count;

    ENTER;
    SAVETMPS;
    PUSHMARK(SP);
    EXTEND(SP, 2);
    PUSHs(pop_data->wrapper);
    PUSHs(my_ga_entity_wrapper(joe));
    PUTBACK;
    count = call_sv(pop_data->evaluate_cb, G_SCALAR);
    SPAGAIN;
    if (count != 1) croak("Callback did not return a value");
    ret_sv = POPs;
    if (SvOK(ret_sv)) {
        ret = 1;
        joe->fitness = SvNV(ret_sv);
    }
    else
        ret = 0;
    PUTBACK;
    FREETMPS;
    LEAVE;
    // fprintf(stderr, "exit cb\n"); fflush(stderr);
    return ret;
}

ga_crossover__perl(population *pop,
                   entity *father, entity *mother,
                   entity *son, entity *daughter) {
    dSP;

    SV *wrapper_father, *wrapper_mother,
       *wrapper_son, *wrapper_daughter;

    if (!father || !mother || !son || !daughter)
        croak("Internal error: null pointer to entity structure passed");

    wrapper_father   = my_ga_entity_wrapper(father);
    wrapper_mother   = my_ga_entity_wrapper(mother);
    wrapper_son      = my_ga_entity_wrapper(son);
    wrapper_daughter = my_ga_entity_wrapper(daughter);

    my_ga_entity_wrapper_set_rw(wrapper_son);
    my_ga_entity_wrapper_set_rw(wrapper_daughter);

    ENTER;
    SAVETMPS;
    PUSHMARK(SP);
    EXTEND(SP, 5);
    PUSHs(pop_data->wrapper);
    PUSHs(wrapper_father);
    PUSHs(wrapper_mother);
    PUSHs(wrapper_son);
    PUSHs(wrapper_daughter);
    PUTBACK;
    call_sv(pop_data->crossover_cb, G_VOID);
    SPAGAIN;
    PUTBACK;
    FREETMPS;
    LEAVE;

    my_ga_entity_wrapper_after_rw(pop, son, wrapper_son);
    my_ga_entity_wrapper_after_rw(pop, daughter, wrapper_daughter);
}

static void
ga__mod__perl(population *pop, entity *joe, SV *perl_cb) {
    dSP;
    SV *wrapper;

    if (!joe) croak("Internal error: null pointer to entity structure passed");

    wrapper  = my_ga_entity_wrapper(joe);
    my_ga_entity_wrapper_set_rw(wrapper);

    ENTER;
    SAVETMPS;
    PUSHMARK(SP);
    EXTEND(SP, 2);
    PUSHs(pop_data->wrapper);
    PUSHs(wrapper);
    PUTBACK;
    call_sv(perl_cb, G_VOID);
    SPAGAIN;
    PUTBACK;
    FREETMPS;
    LEAVE;

    my_ga_entity_wrapper_after_rw(pop, joe, wrapper);
}

static entity *
ga_adapt__perl(population *pop, entity *child) {
    entity *adult = ga_entity_clone(pop, child);
    ga__mod__perl(pop, adult, pop_data->adapt_cb);
    ga_evaluate__perl(pop, adult);
    return adult;
}

static void
ga_mutate__perl(population *pop, entity *father, entity *son) {
    int i;
    if (!father || !son) Perl_croak(aTHX_ "Null pointer to entity structure passed");
    for (i=0; i < pop->num_chromosomes; i++)
        ga_bit_clone(son->chromosome[i], father->chromosome[i], pop->len_chromosomes);
    ga__mod__perl(pop, son, pop_data->mutate_cb);
}

static boolean
ga_seed__perl(population *pop, entity *joe) {
    ga__mod__perl(pop, joe, pop_data->seed_cb);
    return 1;
}

static boolean
ga_stop__perl(int generation, population *pop) {
    dSP;
    int count;
    boolean ret;
    SV *ret_sv;

    ENTER;
    SAVETMPS;
    PUSHMARK(SP);
    EXTEND(SP, 2);
    PUSHs(pop_data->wrapper);
    PUSHs(sv_2mortal(newSVuv(generation)));
    PUTBACK;
    count = call_sv(pop_data->stop_cb, G_SCALAR);
    SPAGAIN;
    if (count != 1)
        croak("Callback did not return a value");
    ret_sv = POPs;
    ret = !SvTRUE(ret_sv);
    PUTBACK;
    FREETMPS;
    LEAVE;
    /*if (!ret) { fprintf(stderr, "stopping\n"); fflush(stderr); } */
    return ret;
}

