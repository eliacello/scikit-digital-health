# Pre-refactor cross-module endpoint comparison

gait_old.gait_endpoints vs gait.gait_metrics on identical per-step inputs.
`identical` = bit-for-bit equal (NaN-aware); `close` = allclose at rtol=atol=1e-10.

| input | old endpoint | new endpoint | delta h | old key | new key | identical | close 1e-10 | max abs diff | note |
|---|---|---|---|---|---|---|---|---|---|
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StrideTime | StrideTime | n/a | `PARAM:stride time` | `stride time` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StrideTime | StrideTime | n/a | `PARAM:stride time asymmetry` | `stride time asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StanceTime | StanceTime | n/a | `PARAM:stance time` | `stance time` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StanceTime | StanceTime | n/a | `PARAM:stance time asymmetry` | `stance time asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | SwingTime | SwingTime | n/a | `PARAM:swing time` | `swing time` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | SwingTime | SwingTime | n/a | `PARAM:swing time asymmetry` | `swing time asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StepTime | StepTime | n/a | `PARAM:step time` | `step time` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StepTime | StepTime | n/a | `PARAM:step time asymmetry` | `step time asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | InitialDoubleSupport | InitialDoubleSupport | n/a | `PARAM:initial double support` | `initial double support` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | InitialDoubleSupport | InitialDoubleSupport | n/a | `PARAM:initial double support asymmetry` | `initial double support asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | TerminalDoubleSupport | TerminalDoubleSupport | n/a | `PARAM:terminal double support` | `terminal double support` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | TerminalDoubleSupport | TerminalDoubleSupport | n/a | `PARAM:terminal double support asymmetry` | `terminal double support asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | DoubleSupport | DoubleSupport | n/a | `PARAM:double support` | `double support` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | DoubleSupport | DoubleSupport | n/a | `PARAM:double support asymmetry` | `double support asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | SingleSupport | SingleSupport | n/a | `PARAM:single support` | `single support` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | SingleSupport | SingleSupport | n/a | `PARAM:single support asymmetry` | `single support asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StepLength | StepLengthModel1 | shared | `PARAM:step length` | `step length m1` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StepLength | StepLengthModel1 | shared | `PARAM:step length asymmetry` | `step length m1 asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StepLength | StepLengthModel1 | native | `PARAM:step length` | `step length m1` | False | False | 1.262e+00 | NaN placement differs |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StepLength | StepLengthModel1 | native | `PARAM:step length asymmetry` | `step length m1 asymmetry` | False | False | 4.463e-01 | NaN placement differs |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StrideLength | StrideLengthModel1 | shared | `PARAM:stride length` | `stride length m1` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StrideLength | StrideLengthModel1 | shared | `PARAM:stride length asymmetry` | `stride length m1 asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StrideLength | StrideLengthModel1 | native | `PARAM:stride length` | `stride length m1` | False | False | 2.450e+00 | NaN placement differs |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StrideLength | StrideLengthModel1 | native | `PARAM:stride length asymmetry` | `stride length m1 asymmetry` | False | False | 3.122e-01 | NaN placement differs |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | GaitSpeed | GaitSpeedModel1 | shared | `PARAM:gait speed` | `gait speed m1` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | GaitSpeed | GaitSpeedModel1 | shared | `PARAM:gait speed asymmetry` | `gait speed m1 asymmetry` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | GaitSpeed | GaitSpeedModel1 | native | `PARAM:gait speed` | `gait speed m1` | False | False | 2.237e+00 | NaN placement differs |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | GaitSpeed | GaitSpeedModel1 | native | `PARAM:gait speed asymmetry` | `gait speed m1 asymmetry` | False | False | 3.154e-01 | NaN placement differs |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | Cadence | Cadence | n/a | `PARAM:cadence` | `cadence` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | IntraStepCovarianceV | IntraStepCovarianceV | n/a | `PARAM:intra-step covariance - V` | `intra-step covariance - V` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | IntraStrideCovarianceV | IntraStrideCovarianceV | n/a | `PARAM:intra-stride covariance - V` | `intra-stride covariance - V` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | HarmonicRatioV | HarmonicRatioV | n/a | `PARAM:harmonic ratio - V` | `harmonic ratio - V` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StrideSPARC | StrideSPARC | n/a | `PARAM:stride SPARC` | `stride SPARC` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | PhaseCoordinationIndex | PhaseCoordinationIndex | n/a | `BOUTPARAM:phase coordination index` | `bout:phase coordination index` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | GaitSymmetryIndex | GaitSymmetryIndex | n/a | `BOUTPARAM:gait symmetry index` | `bout:gait symmetry index` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StepRegularityV | StepRegularityV | n/a | `BOUTPARAM:step regularity - V` | `bout:step regularity - V` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | StrideRegularityV | StrideRegularityV | n/a | `BOUTPARAM:stride regularity - V` | `bout:stride regularity - V` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | AutocovarianceSymmetryV | AutocovarianceSymmetryV | n/a | `BOUTPARAM:autocovariance symmetry - V` | `bout:autocovariance symmetry - V` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | RegularityIndexV | RegularityIndexV | n/a | `BOUTPARAM:regularity index - V` | `bout:regularity index - V` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StrideTime | StrideTime | n/a | `PARAM:stride time` | `stride time` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StrideTime | StrideTime | n/a | `PARAM:stride time asymmetry` | `stride time asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StanceTime | StanceTime | n/a | `PARAM:stance time` | `stance time` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StanceTime | StanceTime | n/a | `PARAM:stance time asymmetry` | `stance time asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | SwingTime | SwingTime | n/a | `PARAM:swing time` | `swing time` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | SwingTime | SwingTime | n/a | `PARAM:swing time asymmetry` | `swing time asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StepTime | StepTime | n/a | `PARAM:step time` | `step time` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StepTime | StepTime | n/a | `PARAM:step time asymmetry` | `step time asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | InitialDoubleSupport | InitialDoubleSupport | n/a | `PARAM:initial double support` | `initial double support` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | InitialDoubleSupport | InitialDoubleSupport | n/a | `PARAM:initial double support asymmetry` | `initial double support asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | TerminalDoubleSupport | TerminalDoubleSupport | n/a | `PARAM:terminal double support` | `terminal double support` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | TerminalDoubleSupport | TerminalDoubleSupport | n/a | `PARAM:terminal double support asymmetry` | `terminal double support asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | DoubleSupport | DoubleSupport | n/a | `PARAM:double support` | `double support` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | DoubleSupport | DoubleSupport | n/a | `PARAM:double support asymmetry` | `double support asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | SingleSupport | SingleSupport | n/a | `PARAM:single support` | `single support` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | SingleSupport | SingleSupport | n/a | `PARAM:single support asymmetry` | `single support asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StepLength | StepLengthModel1 | shared | `PARAM:step length` | `step length m1` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StepLength | StepLengthModel1 | shared | `PARAM:step length asymmetry` | `step length m1 asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StrideLength | StrideLengthModel1 | shared | `PARAM:stride length` | `stride length m1` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StrideLength | StrideLengthModel1 | shared | `PARAM:stride length asymmetry` | `stride length m1 asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | GaitSpeed | GaitSpeedModel1 | shared | `PARAM:gait speed` | `gait speed m1` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | GaitSpeed | GaitSpeedModel1 | shared | `PARAM:gait speed asymmetry` | `gait speed m1 asymmetry` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | Cadence | Cadence | n/a | `PARAM:cadence` | `cadence` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | IntraStepCovarianceV | IntraStepCovarianceV | n/a | `PARAM:intra-step covariance - V` | `intra-step covariance - V` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | IntraStrideCovarianceV | IntraStrideCovarianceV | n/a | `PARAM:intra-stride covariance - V` | `intra-stride covariance - V` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | HarmonicRatioV | HarmonicRatioV | n/a | `PARAM:harmonic ratio - V` | `harmonic ratio - V` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StrideSPARC | StrideSPARC | n/a | `PARAM:stride SPARC` | `stride SPARC` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | PhaseCoordinationIndex | PhaseCoordinationIndex | n/a | `BOUTPARAM:phase coordination index` | `bout:phase coordination index` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | GaitSymmetryIndex | GaitSymmetryIndex | n/a | `BOUTPARAM:gait symmetry index` | `bout:gait symmetry index` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StepRegularityV | StepRegularityV | n/a | `BOUTPARAM:step regularity - V` | `bout:step regularity - V` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | StrideRegularityV | StrideRegularityV | n/a | `BOUTPARAM:stride regularity - V` | `bout:stride regularity - V` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | AutocovarianceSymmetryV | AutocovarianceSymmetryV | n/a | `BOUTPARAM:autocovariance symmetry - V` | `bout:autocovariance symmetry - V` | True | True | 0.000e+00 |  |
| gait pipeline (tests/gait/data/gait_input.npz) | RegularityIndexV | RegularityIndexV | n/a | `BOUTPARAM:regularity index - V` | `bout:regularity index - V` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StrideTime | StrideTime | n/a | `PARAM:stride time` | `stride time` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StrideTime | StrideTime | n/a | `PARAM:stride time asymmetry` | `stride time asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StanceTime | StanceTime | n/a | `PARAM:stance time` | `stance time` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StanceTime | StanceTime | n/a | `PARAM:stance time asymmetry` | `stance time asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | SwingTime | SwingTime | n/a | `PARAM:swing time` | `swing time` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | SwingTime | SwingTime | n/a | `PARAM:swing time asymmetry` | `swing time asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StepTime | StepTime | n/a | `PARAM:step time` | `step time` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StepTime | StepTime | n/a | `PARAM:step time asymmetry` | `step time asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | InitialDoubleSupport | InitialDoubleSupport | n/a | `PARAM:initial double support` | `initial double support` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | InitialDoubleSupport | InitialDoubleSupport | n/a | `PARAM:initial double support asymmetry` | `initial double support asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | TerminalDoubleSupport | TerminalDoubleSupport | n/a | `PARAM:terminal double support` | `terminal double support` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | TerminalDoubleSupport | TerminalDoubleSupport | n/a | `PARAM:terminal double support asymmetry` | `terminal double support asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | DoubleSupport | DoubleSupport | n/a | `PARAM:double support` | `double support` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | DoubleSupport | DoubleSupport | n/a | `PARAM:double support asymmetry` | `double support asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | SingleSupport | SingleSupport | n/a | `PARAM:single support` | `single support` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | SingleSupport | SingleSupport | n/a | `PARAM:single support asymmetry` | `single support asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StepLength | StepLengthModel1 | shared | `PARAM:step length` | `step length m1` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StepLength | StepLengthModel1 | shared | `PARAM:step length asymmetry` | `step length m1 asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StrideLength | StrideLengthModel1 | shared | `PARAM:stride length` | `stride length m1` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StrideLength | StrideLengthModel1 | shared | `PARAM:stride length asymmetry` | `stride length m1 asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | GaitSpeed | GaitSpeedModel1 | shared | `PARAM:gait speed` | `gait speed m1` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | GaitSpeed | GaitSpeedModel1 | shared | `PARAM:gait speed asymmetry` | `gait speed m1 asymmetry` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | Cadence | Cadence | n/a | `PARAM:cadence` | `cadence` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | IntraStepCovarianceV | IntraStepCovarianceV | n/a | `PARAM:intra-step covariance - V` | `intra-step covariance - V` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | IntraStrideCovarianceV | IntraStrideCovarianceV | n/a | `PARAM:intra-stride covariance - V` | `intra-stride covariance - V` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | HarmonicRatioV | HarmonicRatioV | n/a | `PARAM:harmonic ratio - V` | `harmonic ratio - V` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StrideSPARC | StrideSPARC | n/a | `PARAM:stride SPARC` | `stride SPARC` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | PhaseCoordinationIndex | PhaseCoordinationIndex | n/a | `BOUTPARAM:phase coordination index` | `bout:phase coordination index` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | GaitSymmetryIndex | GaitSymmetryIndex | n/a | `BOUTPARAM:gait symmetry index` | `bout:gait symmetry index` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StepRegularityV | StepRegularityV | n/a | `BOUTPARAM:step regularity - V` | `bout:step regularity - V` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | StrideRegularityV | StrideRegularityV | n/a | `BOUTPARAM:stride regularity - V` | `bout:stride regularity - V` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | AutocovarianceSymmetryV | AutocovarianceSymmetryV | n/a | `BOUTPARAM:autocovariance symmetry - V` | `bout:autocovariance symmetry - V` | True | True | 0.000e+00 |  |
| synthetic unit-test fixture | RegularityIndexV | RegularityIndexV | n/a | `BOUTPARAM:regularity index - V` | `bout:regularity index - V` | True | True | 0.000e+00 |  |
| gait_old pipeline (tests/gait_old/data/gait_input.npz) | gait_old.get_strides (internal) | DeltaHModel1 | native | `delta h` | `m1 delta h` | False | False | 2.994e+00 |  |
