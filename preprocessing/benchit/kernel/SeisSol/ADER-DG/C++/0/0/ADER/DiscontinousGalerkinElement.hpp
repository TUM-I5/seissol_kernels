#include "seissol_src/Initializer/preProcessorMacros.fpp"

/**
 * Disconntiunous Galerkin data for a single element.
 */
class DiscontinuousGalerkinElement {
  public:
    //! unknowns $Q_k$
    double unknowns[NUMBEROFUNKNOWNS];
    
    //! star matrices $A^*_k$, $B^*_k$, $C^*_k$
    double aStarFlat[STARMATRIX_NUMBEROFNONZEROS];
    double bStarFlat[STARMATRIX_NUMBEROFNONZEROS];
    double cStarFlat[STARMATRIX_NUMBEROFNONZEROS];
    
    //! time integrated unknowns $\int_{t^n}^{t^{n+1}} Q_k \d t$
    double timeIntegratedUnknowns[NUMBEROFUNKNOWNS];
  
    // Constructor
    DiscontinuousGalerkinElement() {
      // initialize unknowns and time integrated unknowns.
      for( int l_i = 0; l_i < NUMBEROFUNKNOWNS; l_i++) {
        unknowns[l_i] = 1.0;
        timeIntegratedUnknowns[l_i] = 1.0;
      }
      
      // initialize star matrices.
      for( int l_i = 0; l_i < STARMATRIX_NUMBEROFNONZEROS; l_i++) {
        aStarFlat[l_i] = 1.0;
        bStarFlat[l_i] = 1.0;
        cStarFlat[l_i] = 1.0;
      }
    }
};
