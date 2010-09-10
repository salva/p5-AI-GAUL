
static boolean
ga_evaluate__perl(population *pop, entity *joe) {
    dSP;
    int i;
    int count;
    boolean ret;
    SV *ret_sv;
    int bytes = pop_data->chromosome_bytes;
    SV *perl_cb = pop_data->evaluate_cb;
    byte **chromo = (byte **)joe->chromosome;

    if (!chromo) croak("bad entity");
    ENTER;
    SAVETMPS;
    PUSHMARK(SP);
    XPUSHs(sv_2mortal(newSVsv(pop_data->wrapper)));
    for (i = 0; i < pop->num_chromosomes; i++) {
        SV *sv = get_chromosome_sv(pop, i);
        sv_setpvn(sv, chromo[i], bytes);
        XPUSHs(sv);
    }
    PUTBACK;
    count = call_sv(perl_cb, G_SCALAR);
    SPAGAIN;
    if (count != 1) croak("callback did not return a value");
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

static boolean
ga__mod__perl(population *pop, entity *joe, SV *perl_cb) {
    dSP;
    int i;
    byte **chromo = (byte **)joe->chromosome;
    SV **svs = pop_data->chromosome_svs;
    int bytes = pop_data->chromosome_bytes;

    if (!chromo) croak("bad entity");
    ENTER;
    SAVETMPS;
    PUSHMARK(SP);
    XPUSHs(sv_2mortal(newSVsv(pop_data->wrapper)));
    for (i = 0; i < pop->num_chromosomes; i++) {
        SV *sv = get_chromosome_sv(pop, i);
        sv_setpvn(sv, chromo[i], bytes);
        XPUSHs(sv);
    }
    // XPUSHs(sv_2mortal(newSV(0)));
    PUTBACK;
    call_sv(perl_cb, G_SCALAR);
    SPAGAIN;
    PUTBACK;
    FREETMPS;
    LEAVE;
    
    for (i = 0; i < pop->num_chromosomes; i++) {
        STRLEN len;
        const char *pv = SvPV_const(svs[i], len);
        if (len != bytes)
            croak("invalid size %d for adapted chromosome, %d expected", len, bytes);
        Copy(pv, chromo[i], bytes, char);
    }
    return 1;
}

static entity *
ga_adapt__perl(population *pop, entity *child) {
    entity *adult = ga_entity_clone(pop, child);
    ga__mod__perl(pop, adult, pop_data->adapt_cb);
    ga_evaluate__perl(pop, adult);
    return adult;
}

static boolean
ga_seed__perl(population *pop, entity *joe) {
    ga__mod__perl(pop, joe, pop_data->seed_cb);
}

static boolean
ga_stop__perl(int generation, population *pop) {
    dSP;
    int i;
    int count;
    boolean ret;
    SV *ret_sv;
    int bytes = pop_data->chromosome_bytes;
    SV *perl_cb = pop_data->stop_cb;

    ENTER;
    SAVETMPS;
    PUSHMARK(SP);
    XPUSHs(sv_2mortal(newSVsv(pop_data->wrapper)));
    XPUSHs(sv_2mortal(newSVuv(generation)));
    PUTBACK;
    count = call_sv(perl_cb, G_SCALAR);
    SPAGAIN;
    if (count != 1) croak("callback did not return a value");
    ret_sv = POPs;
    // sv_dump(ret_sv);
    ret = !SvTRUE(ret_sv);
    PUTBACK;
    FREETMPS;
    LEAVE;
    // fprintf(stderr, "exit cb\n"); fflush(stderr);
    if (!ret) { fprintf(stderr, "stopping\n"); fflush(stderr); }
    return ret;
}
