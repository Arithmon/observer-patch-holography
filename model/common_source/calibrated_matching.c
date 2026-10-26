/* Forward matching through the separately installed SMDR 1.3 library.
 * This adapter supplies no OPH derivation or native interacting evolution.
 * SMDR: S. P. Martin and D. G. Robertson, arXiv:1907.02500.
 */
#include <stdio.h>
#include "smdr.h"

#define EMIT(name) printf("\"" #name "\":\"%.21Lg\",", (long double) name)

int main(int argc, char **argv) {
    if (argc != 2) return 2;
    char *inputs[] = {
        "SMDR_Q_in", "SMDR_v_in", "SMDR_lambda_in", "SMDR_g3_in",
        "SMDR_g_in", "SMDR_gp_in", "SMDR_yt_in", "SMDR_yb_in",
        "SMDR_yc_in", "SMDR_ys_in", "SMDR_yu_in", "SMDR_yd_in",
        "SMDR_ytau_in", "SMDR_ymu_in", "SMDR_ye_in",
        "SMDR_Delta_alpha_had_5_MZ_in"
    };
    SMDR_Read_Values(argv[1], 16, inputs);
    SMDR_Load_Inputs();
    /* SMDR 1.3 Load_Inputs does not refresh derived coupling powers, and
     * RGrun_SM returns before Update on a zero-length step. Run the common
     * input point to 200 GeV once so none of the downstream matching scales
     * equals the saved input scale. No physical parameter is refitted here.
     */
    SMDR_Update();
    SMDR_RGrun_SM(200., 5);
    SMDR_Save_Inputs();
    SMDR_m2_in = SMDR_m2 = SMDR_Lambda_in = SMDR_Lambda = 0;

    SMDR_Eval_Mt_pole(SMDR_Mt_EXPT, 1, 4, 2, &SMDR_Mt_pole, &SMDR_Gammat_pole);
    SMDR_Eval_Mh_pole(160., 2.5, &SMDR_Mh_pole, &SMDR_Gammah_pole);
    SMDR_Eval_MZ_pole(160., 2.5, &SMDR_MZ_pole, &SMDR_GammaZ_pole,
                      &SMDR_MZ_PDG, &SMDR_GammaZ_PDG);
    SMDR_Eval_MW_pole(160., 2.5, &SMDR_MW_pole, &SMDR_GammaW_pole,
                      &SMDR_MW_PDG, &SMDR_GammaW_PDG);
    SMDR_GFermi = SMDR_Eval_GFermi(SMDR_Mt_EXPT, 3);
    SMDR_Eval_Gauge(SMDR_Mt_pole, SMDR_Mh_pole, SMDR_MW_PDG);
    SMDR_Eval_QCDQED_at_MZ(SMDR_MZ_EXPT, SMDR_MZ_EXPT, 5);
    SMDR_mbmb = SMDR_Eval_mbmb(SMDR_MZ_EXPT, 5);
    SMDR_mcmc = SMDR_Eval_mcmc(SMDR_Mtau_EXPT, SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 5);
    SMDR_Eval_mquarks_2GeV(SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 5,
                         &SMDR_ms_2GeV, &SMDR_mu_2GeV, &SMDR_md_2GeV);
    SMDR_Mtau_pole = SMDR_Eval_Mtau_pole(SMDR_Mtau_EXPT, SMDR_mbmb_EXPT,
                                      SMDR_MZ_EXPT, 3);
    SMDR_Mmuon_pole = SMDR_Eval_Mmuon_pole(SMDR_mcmc_EXPT, SMDR_mcmc_EXPT,
                        SMDR_Mtau_EXPT, SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 2);
    SMDR_Melectron_pole = SMDR_Eval_Melectron_pole(SMDR_mcmc_EXPT, SMDR_mcmc_EXPT,
                        SMDR_Mtau_EXPT, SMDR_mbmb_EXPT, SMDR_MZ_EXPT, 2);

    printf("\nOPH_MATCHING_JSON={");
    EMIT(SMDR_Mt_pole); EMIT(SMDR_Mh_pole);
    EMIT(SMDR_MZ_PDG); EMIT(SMDR_MZ_pole);
    EMIT(SMDR_MW_PDG); EMIT(SMDR_MW_pole);
    EMIT(SMDR_mbmb); EMIT(SMDR_mcmc); EMIT(SMDR_ms_2GeV);
    EMIT(SMDR_mu_2GeV); EMIT(SMDR_md_2GeV);
    EMIT(SMDR_Mtau_pole); EMIT(SMDR_Mmuon_pole); EMIT(SMDR_Melectron_pole);
    EMIT(SMDR_alpha); EMIT(SMDR_alphaS_5_MZ); EMIT(SMDR_GFermi);
    EMIT(SMDR_Gammat_pole); EMIT(SMDR_Gammah_pole);
    EMIT(SMDR_GammaZ_pole); EMIT(SMDR_GammaW_pole);
    printf("\"sizeof_long_double\":%zu}\n", sizeof(long double));
    return 0;
}
