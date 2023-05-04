MODULE ConvEquiResults
    !
    USE globals
    USE generic_routines
    USE QL_routines
    USE EquilibriumCheck
    !
    ! Computes check for best response and equilibrium in all states and for all agents
    !
    IMPLICIT NONE
    !
CONTAINS
    !
    ! &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
    !
    SUBROUTINE ComputeConvEquiResults ( iExperiment )
        !
        ! Computes statistics for one model
        !
        IMPLICIT NONE
        !
        ! Declaring dummy variables
        !
        INTEGER, INTENT(IN) :: iExperiment
        !
        ! Declaring local variable
        !
        INTEGER :: iSession, iAgent, iThres, i, j, iState, &
            OptimalStrategyVec(lengthStrategies), OptimalStrategy(numStates,numAgents), &
            CycleStatesSession(numPeriods), &
            CycleLengthSession, numCycleLength(0:numThresCycleLength), &
            ThresCycleLength(numThresCycleLength)
        INTEGER, DIMENSION(numAgents) :: flagBRAllSession, flagBROnPathSession, flagBROffPathSession
        INTEGER :: flagEQAllSession, flagEQOnPathSession, flagEQOffPathSession
        LOGICAL, DIMENSION(numEquilibriumTypes, numSessions) :: flagEQ
        INTEGER, DIMENSION(numEquilibriumTypes) :: numSessionsByType
        INTEGER :: typeIndex
        !
        REAL(8) :: r_num
        REAL(8), DIMENSION(numAgents) :: freqBRAllSession, freqBROnPathSession, freqBROffPathSession
        REAL(8) :: freqEQAllSession, freqEQOnPathSession, freqEQOffPathSession
        !
        ! Final Profit Results (for all sessions, those that converged to NE, and those that
        ! converged to SPNE)
        !
        REAL(8) :: Profits(numSessions,numAgents), AvgProfits(numSessions)
        REAL(8), DIMENSION(numEquilibriumTypes, numAgents) :: meanProfit, seProfit, meanProfitGain, seProfitGain, &
            NashProfitsByType, CoopProfitsByType
        REAL(8), DIMENSION(numEquilibriumTypes) :: meanAvgProfit, seAvgProfit, meanAvgProfitGain, seAvgProfitGain, &
            MeanNashProfitsByType, MeanCoopProfitsByType
        !
        ! Beginning execution
        !
        PRINT*, 'Computing convergence results given equilibrium results'
        ThresCycleLength = (/ (i, i = 1, numThresCycleLength) /)
        numSessionsByType = 0
        !
        ! Initializing variables
        !
!$      CALL OMP_SET_NUM_THREADS(numCores)
        !
        ! Reading strategies and states at convergence from file
        !
        CALL ReadInfoExperiment()
        !
        ! Beginning loop over sessions
        !
        !$omp parallel do &
        !$omp private(OptimalStrategy,OptimalStrategyVec,CycleLengthSession,CycleStatesSession, &
        !$omp   freqBRAllSession,freqBROnPathSession,freqBROffPathSession,freqEQAllSession,freqEQOnPathSession,freqEQOffPathSession, &
        !$omp   flagBRAllSession,flagBROnPathSession,flagBROffPathSession,flagEQAllSession,flagEQOnPathSession,flagEQOffPathSession) &
        !$omp firstprivate(converged)
        DO iSession = 1, numSessions                  ! Start of loop aver sessions
            !
            PRINT*, 'iSession = ', iSession
            !
            !$omp critical
            OptimalStrategyVec = indexStrategies(:,iSession)
            CycleLengthSession = CycleLength(iSession)
            CycleStatesSession(:CycleLengthSession) = CycleStates(:CycleLengthSession,iSession)
            Profits(iSession,:) = SUM(CycleProfits(:,:CycleLengthSession,iSession),DIM = 2) / &
                DBLE(CycleLengthSession)
            !$omp end critical
            !
            ! Check if the session converged
            !
            IF (converged(iSession) .EQ. 0) THEN
                !$omp critical
                !
                ! This session did not converge. Flag all false
                !
                flagEQ(1, iSession) = .false.
                flagEQ(2, iSession) = .false.
                flagEQ(3, iSession) = .false.
                !
                !$omp end critical
            ELSE
                !
                OptimalStrategy = RESHAPE(OptimalStrategyVec, (/ numStates,numAgents /) )
                !
                CALL computeEqCheckSession(OptimalStrategy,CycleLengthSession,CycleStatesSession(:CycleLengthSession), &
                    freqBRAllSession,freqBROnPathSession,freqBROffPathSession,freqEQAllSession,freqEQOnPathSession,freqEQOffPathSession, &
                    flagBRAllSession,flagBROnPathSession,flagBROffPathSession,flagEQAllSession,flagEQOnPathSession,flagEQOffPathSession)
                !
                !$omp critical
                IF (fixStrategies .EQ. 1) THEN
                    !
                    flagEQ(1, iSession) = .true.
                    flagEQ(2, iSession) = flagBROnPathSession(1) .EQ. 1
                    flagEQ(3, iSession) = flagBRAllSession(1) .EQ. 1
                    numSessionsByType(:) = numSessionsByType(:) + (/ 1, flagBROnPathSession(1), flagBRAllSession(1) /)
                    !
                ELSE 
                    !
                    flagEQ(1, iSession) = .true.
                    flagEQ(2, iSession) = flagEQOnPathSession .EQ. 1
                    flagEQ(3, iSession) = flagEQAllSession .EQ. 1
                    numSessionsByType(:) = numSessionsByType(:) + (/ 1, flagEQOnPathSession, flagEQAllSession /)
                    !
                END IF
                !$omp end critical
                !
            END IF
            !
        END DO                                  ! End of loop over sessions
        !$omp end parallel do
        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Computing averages and descriptive statistics
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        ! Profits
        !
        AvgProfits = SUM(Profits,DIM = 2)/DBLE(numAgents)
        NashProfitsByType = RESHAPE((/ (NashProfits, typeIndex=1,numEquilibriumTypes) /), (/ numEquilibriumTypes, numAgents /))
        CoopProfitsByType = RESHAPE((/ (CoopProfits, typeIndex=1,numEquilibriumTypes) /), (/ numEquilibriumTypes, numAgents /))
        MeanNashProfitsByType = (/ (SUM(NashProfits)/numAgents, typeIndex=1,numEquilibriumTypes)/)
        MeanCoopProfitsByType = (/ (SUM(CoopProfits)/numAgents, typeIndex=1,numEquilibriumTypes)/)

        DO iAgent = 1, numAgents
            !
            DO typeIndex = 1, numEquilibriumTypes
                !
                IF (numSessionsByType(typeIndex) .EQ. 0) THEN
                    !
                    meanProfit(typeIndex, iAgent) = 0.d0
                    seProfit(typeIndex, iAgent) = 0.d0
                    NashProfitsByType(typeIndex, :) = 0.d0
                    !
                ELSE IF (numSessionsByType(typeIndex) .GT. 0) THEN
                    !
                    meanProfit(typeIndex, iAgent) = SUM(Profits(:, iAgent), MASK=flagEQ(typeIndex,:))/DBLE(numSessionsByType(typeIndex))
                    seProfit(typeIndex, iAgent) = SQRT(ABS((SUM(Profits(:,iAgent)**2, MASK=flagEQ(typeIndex,:))/DBLE(numSessionsByType(typeIndex))-meanProfit(typeIndex, iAgent)**2)))
                    !
                END IF
                !
            END DO
            !
        END DO
        DO typeIndex = 1, numEquilibriumTypes
            !
            IF (numSessionsByType(typeIndex) .EQ. 0) THEN
                !
                meanAvgProfit(typeIndex) = 0.d0
                seAvgProfit(typeIndex) = 0.d0
                MeanNashProfitsByType(typeIndex) = 0.d0
                !
            ELSE IF (numSessionsByType(typeIndex) .GT. 0) THEN
                !
                meanAvgProfit(typeIndex) = SUM(AvgProfits(:), MASK=flagEQ(typeIndex,:))/DBLE(numSessionsByType(typeIndex))
                seAvgProfit(typeIndex) = SQRT(ABS((SUM(AvgProfits(:)**2, MASK=flagEQ(typeIndex,:))/DBLE(numSessionsByType(typeIndex))-meanAvgProfit(typeIndex)**2)))
                !
            END IF
            !
        END DO

        meanProfitGain = (meanProfit-NashProfitsByType)/(CoopProfitsByType-NashProfitsByType)
        seProfitGain = seProfit/(CoopProfitsByType-NashProfitsByType)
        meanAvgProfitGain = (meanAvgProfit-MeanNashProfitsByType)/(MeanCoopProfitsByType-MeanNashProfitsByType)
        seAvgProfitGain = seAvgProfit/(MeanCoopProfitsByType-MeanNashProfitsByType)

        !
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Printing averages and descriptive statistics
        ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !
        IF (iExperiment .EQ. 1) THEN
            !
            WRITE(100044,991) &
                ((i, j, i, j, j = 1, numEquilibriumTypes), i = 1, numAgents), &
                (i, i = 1, numEquilibriumTypes), &
                (i, i = 1, numEquilibriumTypes), &
                ((i, j, i, j, j = 1, numEquilibriumTypes), i = 1, numAgents), &
                (i, i = 1, numEquilibriumTypes), &
                (i, i = 1, numEquilibriumTypes), &
                (i, i = 1, numEquilibriumTypes)
991         FORMAT('Experiment ', &
                <numAgents>(<numEquilibriumTypes>('    avgProf', I1, '_', I1, '     seProf', I1, '_', I1)), &
                <numEquilibriumTypes>('     avgProf_', I1), &
                <numEquilibriumTypes>('      seProf_', I1), &
                <numAgents>(<numEquilibriumTypes>('  avgPrGain', I1, '_', I1, '   sePrGain', I1, '_', I1)), &
                <numEquilibriumTypes>('   avgPrGain_', I1), &
                <numEquilibriumTypes>('    sePrGain_', I1), &
                <numEquilibriumTypes>(' numSessions_', I1) &
                )
            !
        END IF
        !
        WRITE(100044,992) codExperiment, &
            ((meanProfit(j, i), seProfit(j, i), j = 1, numEquilibriumTypes), i = 1, numAgents), &
            (meanAvgProfit(i), i = 1, numEquilibriumTypes), &
            (seAvgProfit(i), i = 1, numEquilibriumTypes), &
            ((meanProfitGain(j, i), seProfitGain(j, i), j = 1, numEquilibriumTypes), i = 1, numAgents), &
            (meanAvgProfitGain(i), i = 1, numEquilibriumTypes), &
            (seAvgProfitGain(i), i = 1, numEquilibriumTypes), &
            (numSessionsByType(i), i = 1, numEquilibriumTypes)
992     FORMAT(I10, 1X, &
            <6*(numAgents+1)>(F14.5), &
            <6*(numAgents+1)>(F14.5), &
            (I14), (I14), (I14) &
            )
        !
        ! Ending execution and returning control (numEquilibriumTypes, numAgents)
        !
    END SUBROUTINE ComputeConvEquiResults
END MODULE ConvEquiResults

