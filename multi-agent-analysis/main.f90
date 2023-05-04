PROGRAM main
!
    USE globals
    USE LearningSimulation
    USE ConvergenceResults
    USE ConvEquiResults
    USE ImpulseResponse
    USE EquilibriumCheck
    USE QGapToMaximum
    USE LearningTrajectory
    USE DetailedAnalysis
    USE QL_routines
    USE PI_routines
    USE generic_routines
!
    IMPLICIT NONE
!
! Declaring variables and parameters
!
    INTEGER :: iExperiment, i, iAgent
    CHARACTER(len = 50) :: MainFolder
    INTEGER :: argI
!
! Beginning execution
!
!
! Get command line arguments
!
CALL get_command_argument(1, MainFolder)
IF (LEN_TRIM(MainFolder) == 0) THEN
    PRINT*, "Main folder not specified"
    ERROR STOP
END IF
!
! Opening files
!
    OPEN(UNIT = 10001,FILE = TRIM(MainFolder)//"/A_InputParameters.txt")
    CALL readBatchVariables(10001)
!
    OPEN(UNIT = 10002,FILE = TRIM(MainFolder)//"/A_res.txt")
    OPEN(UNIT = 100022,FILE = TRIM(MainFolder)//"/A_convResults.txt")
    IF (SwitchProfitAnalysisByAction .EQ. 1) OPEN(UNIT = 100033,FILE = TRIM(MainFolder)//"/A_m_results.txt")
    IF (SwitchProfitAnalysisByActionByEquiType .EQ. 1) OPEN(UNIT = 100044,FILE = TRIM(MainFolder)//"/A_m_eq_results.txt")
    IF (SwitchImpulseResponseToBR .EQ. 1) OPEN(UNIT = 10003,FILE = TRIM(MainFolder)//"/A_irToBR.txt")
    IF (SwitchImpulseResponseToNash .GE. 1) OPEN(UNIT = 100031,FILE = TRIM(MainFolder)//"/A_irToNash.txt")
    IF (SwitchImpulseResponseToAll .EQ. 1) OPEN(UNIT = 100032,FILE = TRIM(MainFolder)//"/A_irToAll.txt")
    IF (SwitchEquilibriumCheck .EQ. 1) OPEN(UNIT = 10004,FILE = TRIM(MainFolder)//"/A_ec.txt")
    IF (SwitchQGapToMaximum .EQ. 1) OPEN(UNIT = 10006,FILE = TRIM(MainFolder)//"/A_qg.txt")
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Loop over models
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!
    DO iExperiment = 1, numExperiments
        !
        ! Read model parameters
        !
        CALL readExperimentVariables(10001)
        labelStates = computeStatesCodePrint()
        !
        ! Precompute profit matrix for every action pair
        !
        CALL computePIMatricesLogit(DemandParameters,NashPrices,CoopPrices,&
            PI,NashProfits,CoopProfits,NashMarketShares,CoopMarketShares,PricesGrids)
        !
        ! Precompute profit gains (from nash equilibrium profit) for every action pair
        !
        DO iAgent = 1, numAgents
            !
            PG(:,iAgent) = (PI(:,iAgent)-NashProfits(iAgent))/(CoopProfits(iAgent)-NashProfits(iAgent))
            !
        END DO
        !
        ! Creating I/O filenames
        !
        WRITE(ExperimentNumber, "(I0.<LengthFormatTotExperimentsPrint>, A4)") codExperiment, ".txt"
        FileNameInfoExperiment = TRIM(MainFolder)//"/InfoExperiment_"//ExperimentNumber
        !
        ! Print message
        !
        WRITE(*,11) iExperiment, numExperiments, numCores
11      FORMAT('model = ', I6, ' / numExperiments = ', I6, ' / numCores = ', I6)
        IF (SwitchComputeSimulation .EQ. 1) THEN
            !
            ! Compute QL strategy
            !
            CALL computeExperiment(iExperiment,codExperiment,alpha,ExplorationParameters,delta)
            !
            ! Results at convergence
            !
            CALL ComputeConvResults(iExperiment)
            !
        END IF        
        !
        ! Analyze how profits differ by actions segmented by equilibrium concepts
        !
        IF (SwitchProfitAnalysisByActionByEquiType .EQ. 1) CALL ComputeConvEquiResults(iExperiment)
        !
        ! Impulse Response analysis to one-period deviation to static best response
        ! NB: The last argument in computeIRAnalysis is "IRType", and it's crucial:
        ! IRType < 0 : One-period deviation to the price IRType
        ! IRType = 0 : One-period deviation to static BR
        ! IRType > 0 : IRType-period deviation to Nash
        !
        IF (SwitchImpulseResponseToBR .EQ. 1) CALL computeIRAnalysis(iExperiment,10003,0)
        !
        ! Impulse Response to a permanent or transitory deviation to Nash prices
        !
        IF (SwitchImpulseResponseToNash .GE. 1) CALL computeIRAnalysis(iExperiment,100031,SwitchImpulseResponseToNash)
        !
        ! Impulse Response analysis to one-period deviation to all prices
        !
        IF (SwitchImpulseResponseToAll .EQ. 1) THEN
            !
            DO i = 1, numPrices
                !
                CALL computeIRAnalysis(iExperiment,100032,-i)
                !
            END DO
            !
        END IF
        !
        ! Equilibrium Check
        !
        IF (SwitchEquilibriumCheck .EQ. 1) CALL computeEqCheck(iExperiment)
        !
        ! Q Gap w.r.t. Maximum
        !
        IF (SwitchQGapToMaximum .EQ. 1) CALL computeQGapToMax(iExperiment)
        !
        ! Learning Trajectory analysis
        !
        IF (ParamsLearningTrajectory(1) .GT. 0) &
            CALL ComputeLearningTrajectory(iExperiment,codExperiment,alpha,ExplorationParameters,delta)
        !
        ! Detailed Impulse Response analysis to one-period deviation to all prices
        !
        IF (SwitchDetailedAnalysis .EQ. 1) CALL ComputeDetailedAnalysis(iExperiment)
        !
        ! Deallocating arrays
        !
        CALL closeSession()
        !
        ! End of loop over models
        !
    END DO
!
! Deallocating arrays
!
    CALL closeBatch()
!
! Closing output files
!
    CLOSE(UNIT = 10001)
    CLOSE(UNIT = 10002)
    CLOSE(UNIT = 100022)
    IF (SwitchImpulseResponseToBR .EQ. 1) CLOSE(UNIT = 10003)
    IF (SwitchImpulseResponseToNash .GE. 1) CLOSE(UNIT = 100031)
    IF (SwitchImpulseResponseToAll .EQ. 1) CLOSE(UNIT = 100032)
    IF (SwitchEquilibriumCheck .EQ. 1) CLOSE(UNIT = 10004)
    IF (SwitchQGapToMaximum .EQ. 1) CLOSE(UNIT = 10006)
!
! End of execution
!
END PROGRAM main
