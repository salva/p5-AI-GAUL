
static SV *
pop_av_slot_(population *pop, int slot) {
    SV **sv = av_fetch((AV*)SvRV((SV*)pop->data), slot, 0);
    if (!sv) croak("internal error, callback for slot %d is empty", slot);
    return *sv;
}

#define pop_av_slot(slot) pop_av_slot_(pop, slot)

static boolean ga_evaluate__perl(population *pop, entity *joe) {
    dSP;
    int count;
    boolean ret;
    SV *ret_sv;
    // fprintf(stderr, "enter cb\n"); fflush(stderr);
    SV *perl_cb = pop_data->evaluate_cb;
    ENTER;
    SAVETMPS;
    PUSHMARK(SP);
    XPUSHs(sv_2mortal(newSVsv(pop->data)));
    XPUSHs(sv_2mortal(entity_to_sv(pop, joe)));
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
