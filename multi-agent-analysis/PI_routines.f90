MODULE PI_routines
!
    USE globals
    USE generic_routines
    USE QL_routines
!
! Various routines used to compute PI matrices at runtime
!
    IMPLICIT NONE
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE computePIMatricesLogit ( DemandParameters, NashPrices, CoopPrices, &
        PI, NashProfits, CoopProfits, &
        NashMarketShares, CoopMarketShares, &
        PricesGrids )
        !
        ! Computes the Logit common payoff matrix PI
        !
        IMPLICIT NONE
        !
        ! Declaring dummy variables
        !
        REAL(8), INTENT(IN) :: DemandParameters(numDemandParameters)
        REAL(8), DIMENSION(numAgents), INTENT(IN) :: NashPrices, CoopPrices
        REAL(8), INTENT(OUT) :: PI(numActions,numAgents)
        REAL(8), DIMENSION(numAgents), INTENT(OUT) :: NashProfits, CoopProfits, &
            NashMarketShares, CoopMarketShares
        REAL(8), DIMENSION(numPrices,numAgents), INTENT(OUT) :: PricesGrids
        !
        ! Declaring local variables
        !
        REAL(8) :: a0, mu, extend(2)
        REAL(8), DIMENSION(numAgents) :: a, c, d, stepPrices, prices
        INTEGER :: i, j, iter, iAgent
        !
        ! Beginning execution
        !
        ! Extract demand parameters
        !
        a0 = DemandParameters(1)
        a = DemandParameters(2:1+numAgents)
        c = DemandParameters(2+numAgents:1+2*numAgents)
        mu = DemandParameters(2+2*numAgents)
        extend = DemandParameters(3+2*numAgents:4+2*numAgents)
        !
        ! 1. Compute repeated Nash profits
        !
        NashMarketShares = logitDemands(a0,a,c,mu,NashPrices)
        NashProfits = (NashPrices-c)*NashMarketShares
        !
        ! 2. Compute cooperation profits
        !
        CoopMarketShares = logitDemands(a0,a,c,mu,CoopPrices)
        CoopProfits = (CoopPrices-c)*CoopMarketShares
        !
        ! 3. Compute price grid
        !
        ! Upper and lower bounds
        !
        IF (ALL(extend .GT. 0.d0)) THEN
            !
            ! Lower bound = pNash - extend(1)*(pCoop - pNash)
            ! Upper bound = pCoop + extend(2)*(pCoop - pNash)
            !
            PricesGrids(1,:) = NashPrices-extend(1)*(CoopPrices-NashPrices)
            PricesGrids(numPrices,:) = CoopPrices+extend(2)*(CoopPrices-NashPrices)
            !
        ELSE IF ((extend(1) .LT. 0.d0) .AND. (extend(2) .GE. -EPSILON(extend(2)))) THEN
            !
            ! Lower bound = cost + extend(1)*cost
            ! Upper bound = pCoop + extend(2)*(CoopPrices-NashPrices)
            !
            PricesGrids(1,:) = c+extend(1)*c
            PricesGrids(numPrices,:) = CoopPrices+extend(2)*(CoopPrices-NashPrices)
            !
        END IF
        !
        ! Grids
        !
        stepPrices = (PricesGrids(numPrices,:)-PricesGrids(1,:))/(numPrices-1)
        DO i = 2, numPrices-1
            !
            PricesGrids(i,:) = PricesGrids(i-1,:)+stepPrices
            !
        END DO
        !
        ! 4. Compute Pi matrices
        !
        DO i = 1, numActions
            !
            DO j = 1, numAgents
                !
                prices(j) = PricesGrids(indexActions(i,j),j)
                !
            END DO
            !
            d = logitDemands(a0,a,c,mu,prices)
            PI(i,:) = (prices-c)*d
            !
        END DO
        !
        ! Ending execution and returning control
        !
    END SUBROUTINE computePIMatricesLogit
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    FUNCTION logitDemands ( a0, a, c, mu, p )
        !
        ! Computes logit demands
        !
        IMPLICIT NONE
        !
        ! Declaring dummy variables
        !
        REAL(8), INTENT(IN) :: a0, mu
        REAL(8), DIMENSION(numAgents), INTENT(IN) :: a, c, p
        !
        ! Declaring function's type
        !
        REAL(8), DIMENSION(numAgents) :: logitDemands
        !
        ! Beginning execution
        !
        logitDemands = EXP((a-p)/mu)
        logitDemands = logitDemands/(SUM(logitDemands)+EXP(a0/mu))
        !
        ! Ending execution and returning control
        !
    END FUNCTION logitDemands
!
! End of execution
!
END MODULE PI_routines
