#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include "siena07setup.h"
#include "siena07models.h"
#include <R_ext/Rdynload.h>

//using namespace std;

/*
  From tools::package_native_routine_registration_skeleton.
  And I looked at how it is done in base package methods.
*/

/* .Call calls */

// I tried various definitions of the arguments, but none worked:
extern SEXP C_Behavior(SEXP RpData, SEXP BEHLIST);
//extern SEXP C_Bipartite(SEXP, SEXP);
extern SEXP C_Bipartite(void *);
//extern SEXP C_ChangingCovariates(SEXP, SEXP);
extern SEXP C_ChangingCovariates();
extern SEXP C_ChangingDyadicCovariates(SEXP, SEXP);
extern SEXP C_ConstantCovariates(SEXP, SEXP);
extern SEXP C_Constraints(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_deleteData(SEXP);
extern SEXP C_deleteModel(SEXP);
extern SEXP C_DyadicCovariates(SEXP, SEXP);
extern SEXP C_effects(SEXP, SEXP);
extern SEXP C_ExogEvent(SEXP, SEXP);
extern SEXP C_forwardModel(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_getTargets(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_interactionEffects(SEXP, SEXP);
extern SEXP C_mlInitializeSubProcesses(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_mlMakeChains(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_mlPeriod(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP C_OneMode(SEXP, SEXP);
extern SEXP C_setupData(SEXP, SEXP);
extern SEXP C_setupModelOptions(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);


#define CALLDEF(name, n)  {#name, (DL_FUNC) &name, n}

static const R_CallMethodDef CallEntries[] = {
   CALLDEF(Behavior, 2),
   CALLDEF(Bipartite, 2),
   CALLDEF(ChangingCovariates, 2),
   CALLDEF(ChangingDyadicCovariates, 2),
   CALLDEF(ConstantCovariates, 2),
   CALLDEF(Constraints, 7),
   CALLDEF(deleteData, 1),
   CALLDEF(deleteModel, 1),
   CALLDEF(DyadicCovariates, 2),
   CALLDEF(effects, 2),
   CALLDEF(ExogEvent, 2),
   CALLDEF(forwardModel, 16),
   CALLDEF(getTargets, 6),
   CALLDEF(interactionEffects, 2),
   CALLDEF(mlInitializeSubProcesses, 10),
   CALLDEF(mlMakeChains, 9),
   CALLDEF(mlPeriod, 13),
   CALLDEF(OneMode, 2),
   CALLDEF(setupData, 2),
   CALLDEF(setupModelOptions, 12),
    {NULL, NULL, 0}
};

void R_init_RSiena(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
