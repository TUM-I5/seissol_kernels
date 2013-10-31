#include <iostream>

#include "DiscontinousGalerkinElement.hpp"

/**
 * Benchmark data necessary for CauchyKovalewski / ADER benchmark kernels.
 */
class BenchmarkData {
  public:
    //! minimum  number of elemnts in the pseudo-discetization
    int minimumNumberOfElements;
    
    //! maximum  number of elemnts in the pseudo-discetization
    int maximumNumberOfElements;
    
    //! increment in the number of elements
    int incrementPerStep;
    
    //! total number of steps (different elements numbers in the pseudo discretization)
    int numberOfSteps;
    
    //! transposed stiffness matrices (divided by the mass matrix) $K_\xi^T * M^{-1}$, $K_\eta^T * M^{-1}$, $K_\zeta^T * M^{-1}$
    double kXiDivMTFlat  [KXI_NUMBEROFNONZEROS  ];
    double kEtaDivMTFlat [KETA_NUMBEROFNONZEROS ];
    double kZetaDivMTFlat[KZETA_NUMBEROFNONZEROS];
    
    //! time step width $\Delta t$
    double deltaT;
  
    //! pointer to the discontinuous galerkin elements, which hold the element local data: $Q_k$, $A^*_k$, $B^*_k$, $C^*_k$, $\int_{t^n}^{t^{n+1}} Q_k \d t$
    DiscontinuousGalerkinElement* discontinousGalerkinElements;
    
    //! bool, which keeps track if the element local memory has been alocated
    bool memoryAllocated;
    
    // Constructor
    BenchmarkData() {
      // initialize benchmark information with invalid values.
      minimumNumberOfElements = -1;
      maximumNumberOfElements = -1;
      incrementPerStep        = -1;
      numberOfSteps           = -1;
      
      // initialize time step width.
      deltaT = 1.0;
      
      // initialize stiffness matrices with some values
      for( int l_i = 0; l_i < KXI_NUMBEROFNONZEROS; l_i++  ) {
        kXiDivMTFlat[l_i] = 1.0;
      }
      for( int l_i = 0; l_i < KETA_NUMBEROFNONZEROS; l_i++ ) {
        kEtaDivMTFlat[l_i] = 0.1;
      }
      for( int l_i = 0; l_i < KZETA_NUMBEROFNONZEROS; l_i++) {
        kZetaDivMTFlat[l_i] = 0.3;
      }
      
      // no memory allocated so far
      memoryAllocated = false;
    }
    
    // destructor ( frees element local memory )
    ~BenchmarkData() {
      if( memoryAllocated) {
        delete[] discontinousGalerkinElements;
      }
    }
    
    /**
     * Initialize cell local data of the pseudo-discritization
     * Remark: The same grid is used updates with less than the maximum number of elements appearing in the overall benchmark.
     *
     * @param i_numberOfElements number of elements in the pseudo-discretization.
     */
    void initializeDiscontinousGalerkinElements( int i_numberOfElements ) {
      //std::cout << std::endl << "BenchmarkData: Allocating " << i_numberOfElements << " elements" << std::endl;
      
      // allocate memory for DG elements
      discontinousGalerkinElements = new DiscontinuousGalerkinElement[i_numberOfElements];
      
      // memory has been allocated
      memoryAllocated = true;
    }
};
