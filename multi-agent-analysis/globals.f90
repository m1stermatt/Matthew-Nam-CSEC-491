MODULE globals
!
    USE generic_routines
!
! Declare global parameters and variables
!
    IMPLICIT NONE
!
! Parameters
!
    INTEGER, PARAMETER :: numShockPeriodsPrint = 10
    INTEGER, PARAMETER :: numThresCycleLength = 10
    INTEGER, PARAMETER :: numEquilibriumTypes = 3
!
! Variables
!
    INTEGER :: numExperiments, totExperiments, numCores, numSessions, itersPerEpisode, maxNumEpisodes, maxIters, &
        itersInPerfMeasPeriod, printQ, codExperiment, numPrices, decayLearning, decayExploration, fixStrategies, &
        DepthState0, DepthState, LengthStates, numStates, lengthStrategies, &
        LengthFormatStatesPrint, LengthFormatActionPrint, LengthFormatTotExperimentsPrint,LengthFormatNumSessionsPrint, &
        numAgents, numActions, numDemandParameters, numPeriods, &
        numExplorationParameters, &
        SwitchComputeSimulation, &
        SwitchProfitAnalysisByAction, SwitchProfitAnalysisByActionByEquiType, &
        SwitchImpulseResponseToBR, SwitchImpulseResponseToNash, SwitchImpulseResponseToAll, &
        SwitchEquilibriumCheck, SwitchQGapToMaximum, ParamsLearningTrajectory(2), &
        SwitchDetailedAnalysis
    REAL(8) :: delta, PerfMeasPeriodLength, meanNashProfit, meanCoopProfit, gammaSinghVives
    CHARACTER(len = 100) :: ExperimentNumber, FileNameInfoExperiment
!
    INTEGER, ALLOCATABLE :: converged(:), indexStrategies(:,:), indexLastState(:,:), CycleLength(:), &
        CycleStates(:,:), CyclePrices(:,:,:), &
        indexActions(:,:), cStates(:), cActions(:)
    REAL(8), ALLOCATABLE :: timeToConvergence(:), CycleProfits(:,:,:), &
        NashProfits(:), CoopProfits(:), maxValQ(:,:), NashPrices(:), CoopPrices(:), &
        PI(:,:), PG(:,:), alpha(:), DiscountFactors(:), &
        DemandParameters(:), MExpl(:), ExplorationParameters(:), &
        NashMarketShares(:), CoopMarketShares(:), PricesGrids(:,:)
    CHARACTER(len = :), ALLOCATABLE :: labelStates(:)
    CHARACTER(len = :), ALLOCATABLE :: QFileFolderName(:)
!
CONTAINS
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE readBatchVariables ( unitNumber )
        !
        ! Reads input variables
        !
        IMPLICIT NONE
        !
        ! Declaring dummy variables
        !
        INTEGER, INTENT(IN) :: unitNumber
        !
        ! Beginning execution
        !
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) numExperiments, totExperiments
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) numCores
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) numSessions
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) itersPerEpisode
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) maxNumEpisodes
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) PerfMeasPeriodLength
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) numAgents
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) DepthState0
        DepthState = MAX(1,DepthState0)          ! Accomodates the DepthState = 0 case
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) decayLearning
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) decayExploration
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) fixStrategies
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) SwitchComputeSimulation
        
        !
        ! Global variables
        !
        LengthFormatTotExperimentsPrint = 1+INT(LOG10(DBLE(totExperiments)))
        LengthFormatNumSessionsPrint = 1+INT(LOG10(DBLE(numSessions)))
        maxIters = maxNumEpisodes*itersPerEpisode
        itersInPerfMeasPeriod = INT(PerfMeasPeriodLength*itersPerEpisode)
        LengthStates = MAX(1,numAgents*DepthState0)
        numExplorationParameters = numAgents
        numDemandParameters = 2*numAgents+4 ! a0, ai, ci, mu, extend

        READ(unitNumber,'(1X)')
        !
        ! Continue reading input settings
        !
        READ(unitNumber,*) SwitchProfitAnalysisByAction
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) SwitchProfitAnalysisByActionByEquiType
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) SwitchImpulseResponseToBR
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) SwitchImpulseResponseToNash
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) SwitchImpulseResponseToAll
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) SwitchEquilibriumCheck
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) SwitchQGapToMaximum
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) ParamsLearningTrajectory
        READ(unitNumber,'(1X)')
        READ(unitNumber,*) SwitchDetailedAnalysis
        READ(unitNumber,'(1X)')
        !
        ! Allocating matrices and vectors
        !
        ALLOCATE(converged(numSessions),timeToConvergence(numSessions),indexLastState(LengthStates,numSessions), &
            CycleLength(numSessions), DemandParameters(numDemandParameters), &
            ExplorationParameters(numExplorationParameters), MExpl(numExplorationParameters), &
            cStates(LengthStates),cActions(numAgents), &
            alpha(numAgents),NashProfits(numAgents),CoopProfits(numAgents), &
            NashPrices(numAgents),CoopPrices(numAgents), &
            NashMarketShares(numAgents),CoopMarketShares(numAgents))
        ALLOCATE(CHARACTER(len = 200) :: QFileFolderName(numAgents))
        !
        ! Ending execution and returning control
        !
    END SUBROUTINE readBatchVariables
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE closeBatch ( )
        !
        ! Reads input variables
        !
        IMPLICIT NONE
        !
        ! Beginning execution
        !
        DEALLOCATE(converged,timeToConvergence,indexLastState, &
            CycleLength,NashProfits,CoopProfits,QFileFolderName, &
            alpha,MExpl,ExplorationParameters,NashPrices,CoopPrices, &
            DemandParameters,cStates,cActions, &
            NashMarketShares,CoopMarketShares)
        !
        ! Ending execution and returning control
        !
    END SUBROUTINE closeBatch
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE readExperimentVariables ( unitNumber )
        !
        ! Reads input variables
        !
        IMPLICIT NONE
        !
        ! Declaring dummy variables
        !
        INTEGER, INTENT(IN) :: unitNumber
        !
        ! Declaring local variables
        !
        INTEGER :: i, iAction
        !
        ! Beginning execution
        !
        READ(unitNumber,*) codExperiment, printQ, alpha, MExpl, delta, &
            DemandParameters, NashPrices, CoopPrices, numPrices

        !
        ! Global variables
        !
        LengthFormatStatesPrint = LengthStates*(1+FLOOR(LOG10(DBLE(numPrices))))+LengthStates-1
        numStates = numPrices**(numAgents*DepthState0)
        numPeriods = numStates+1
        numActions = numPrices**numAgents   ! Actions contain combinations of prices;
        lengthStrategies = numAgents*numStates
        lengthFormatActionPrint = FLOOR(LOG10(DBLE(numPrices)))+1
        !
        ! Allocating matrices and vectors
        !
        ALLOCATE(indexStrategies(lengthStrategies,numSessions), &
            CycleStates(numPeriods,numSessions),CyclePrices(numAgents,numPeriods,numSessions),CycleProfits(numAgents,numPeriods,numSessions), &
            indexActions(numActions,numAgents), &
            DiscountFactors(0:numStates), &
            maxValQ(numStates,numAgents), &
            PI(numActions,numAgents),PG(numActions,numAgents), &
            PricesGrids(numPrices,numAgents))
        ALLOCATE(CHARACTER(len = 3+LengthFormatStatesPrint) :: labelStates(numStates))
        !
        cStates = (/ (numPrices**i, i = LengthStates-1, 0, -1) /)
        cActions = (/ (numPrices**i, i = numAgents-1, 0, -1) /)
        !
        ! Actions contain the most recent prices of all agents. Actions and States
        ! coincide when DepthState == 1
        !
        DO iAction = 1, numActions
            !
            indexActions(iAction,:) = convertNumberBase(iAction-1,numPrices,numAgents)
            !
        END DO
        !
        ! Set exploration parameters at scale of iterations
        !
        IF (decayExploration .EQ. 1) THEN
            !
            ExplorationParameters = EXP(-MExpl/DBLE(itersPerEpisode))
            !
        ELSE
            !
            ExplorationParameters = MExpl
            !
        END IF
        !
        ! Sanity check: with memoryless agents, only delta = 0 makes sense
        ! If necessary, any value of delta provided on input is overridden
        !
        IF (DepthState0 .EQ. 0) delta = 0.d0
        !
        ! Compute array of discount factors
        !
        DiscountFactors = (/ (delta**i, i = 0, numPeriods-1) /)
        !
        ! Ending execution and returning control
        !
    END SUBROUTINE readExperimentVariables
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
    SUBROUTINE closeSession ( )
        !
        ! Reads input variables
        !
        IMPLICIT NONE
        !
        ! Beginning execution
        !
        DEALLOCATE(indexStrategies,indexActions,CycleStates,CyclePrices,CycleProfits, &
            DiscountFactors, maxValQ, PI,PG, PricesGrids, labelStates)
        !
        ! Ending execution and returning control
        !
    END SUBROUTINE closeSession
!
! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
END MODULE globals
