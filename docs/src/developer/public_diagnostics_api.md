# Публичные наблюдения и отдельные production operators

Эти in-memory entry points принадлежат Core. Они позволяют Runner читать диагностику и проверять
операторы без зависимости от private modules. Файлы, формат, run/execution identity и привязка к
manifest остаются у Runner/Contracts. Новые записи не входят в persisted exact-state graph. Политики
[точного восстановления](public_state_api.md) и существующий `solution_scientific_assessment`
сохраняются.

## Валидация входов и оценка памяти

```julia
validate_options(x::SolverOptions)::SolverOptions
validate_options(x::AlgorithmOptions)::AlgorithmOptions
validate_options(x::ProductionOptions)::ProductionOptions
validate_options(x::ProductionKernelOptions)::ProductionKernelOptions
validate_result_classification(converged::Bool, status::Symbol, scba_quality::Symbol)::Symbol
estimate_production_memory(n::NumericalParameters, scattering::ScatteringOptions;
    options::ProductionOptions=ProductionOptions(),
    worker_capacity::Union{Nothing,Integer}=nothing)::ProductionMemoryEstimate
```

`validate_options` возвращает тот же объект либо прежний `ArgumentError` owning validator. Не строит
grids/basis/kernels, не вызывает callbacks. Проверка классификации проверяет согласованность трёх
labels и возвращает quality: strict inner SCBA может сочетаться с unconverged outer loop. Это не
измерение истинности labels и не scientific acceptance.

Новый overload памяти использует прежнюю conservative prebuild classification: все enabled
mechanisms входят в count; LO/impurity/IFR считаются dense, acoustic/alloy reduced. После build
остаётся существующая оценка фактических kernels. Overload только оценивает и может вернуть peak
выше budget. Вызывающий Runner обязан проверить resolved budget до build; Core production preflight
по-прежнему выполняет ту же проверку с прежним сообщением ошибки. RSS и runtime этим не измеряются.

## Metadata и positivity

```julia
solution_metric_diagnostics(solution::NEGFSolution)::SolutionMetricDiagnostics
positivity_diagnostics(solution::NEGFSolution;
    integrate::Bool=false)::PositivityDiagnostics
```

`SolutionMetricDiagnostics` содержит detached `summary_metrics`, `strict_limits` и
`approximate_limits::Union{Nothing,Dict{Symbol,Float64}}`. Final report имеет приоритет над history;
dynamic collision/candidate metric names получают прежние limits. Approximate limits появляются
только при enabled diagnostic policy и описывают reporting band. Они не меняют original receipt,
eligibility, convergence или strict acceptance. Getter не запускает fresh map, SCBA или Poisson.
Missing metric отсутствует, NaN/Inf сохраняются. Наличие threshold не означает измеренный pass;
measurement/applicability/scientific flags следует читать из `solution_scientific_assessment`.

`PositivityDiagnostics` содержит `metric_version`, `matrix_kinds`, шесть immutable `witnesses`,
`strict_ratio_limit`, шесть `strict_absolute_limits_scaled`, `strict_threshold_satisfied`,
`scan_order` и optional `integrals`. Порядок kinds фиксирован: `spectral`, `occupied`, `unoccupied`,
`broadening`, `raw_occupied`, `raw_unoccupied`. Последние две матрицы вычисляются как
`GR*(-iΣ<)*GR'` и `GR*(iΣ>)*GR'`; нормировка charge их не исправляет. Scattering суммируется в
owning mechanism order.

Каждый `PositivityWitness` сохраняет status, matrix kind, E/k indices, minimum eigenvalue, block
norm, floor, ratio, absolute/relative defect, backward error, Hermiticity defect и construction
scale. `status=:measured` допускает численное сравнение; `:error` исключает pass. В частности,
private compatibility constructor при nonfinite input заполняет `hermiticity_defect=0`; публичный
результат сохраняет raw значение, но это **не измеренная нулевая defect**. NaN/Inf не заменяются.
Валидатор shapes отсекает несовместимые оси до обхода, failed scientific values остаются
наблюдаемыми.

Light mode сканирует E→k, `integrals=nothing` означает «не запрошено». Full mode сканирует k→E и
переиспользует witnesses этого же прохода. Оба сохраняют strict `>` selection: первый из равных
ratios остаётся witness. Поэтому `scan_order` равен соответственно `:energy_momentum` или
`:momentum_energy`; равные maxima могут иметь разные coordinates. Это прежнее поведение light Core
scan и полного Runner matrix audit.

Full `MatrixIntegralDiagnostics` возвращает следующие owned arrays и immutable counts:

| Поле                                 | Оси                                  | Нормировка / единицы                                                |
| ------------------------------------ | ------------------------------------ | ------------------------------------------------------------------- |
| `spectral_integral`                  | k,a,b                                | Весь complex integral `Σe wE[e] A[e,k]/(2π)`, dimensionless         |
| `spectral_integral_eigenvalues`      | k,state                              | Eigenvalues полного integral, dimensionless                         |
| `spectral_trace_cumulative`          | E,k                                  | Real trace prefix с full-grid weights, dimensionless                |
| `negative_energy_integral_scaled`    | kind,k                               | `Σe wE[e] Σj max(0,-λj)/(2π)`; solver units, broadening — energy²   |
| `negative_block_count`               | kind, NTuple{6,Int}                  | Число блоков с minimum eigenvalue < 0                               |
| `beyond_backward_error_block_count`  | kind, NTuple{6,Int}                  | Minimum eigenvalue < −budget                                        |
| `invalid_block_count`                | kind, NTuple{6,Int}                  | Nonfinite/eigensolve failure; соответствующий negative integral NaN |
| `total_blocks_per_kind`              | scalar                               | NE × Nk                                                             |
| `negative_phase_space_weight_per_m2` | kinds (1,2,3,5,6), NTuple{5,Float64} | `gs Σk wk[k] negative[k]/L0²`, m⁻²                                  |

Broadening не получает density units. Negative eigenvalues не clipping и не вычитаются из state. PSD
metric version 3 и budget `64*n*eps(Float64)*max(norm,construction_scale)+construction_error` не
меняются. Light/full локальный флаг равен measured && `ratio < r_PSD`; display absolute limit равен
`backward_error + r_PSD*block_norm`. Это не global acceptance. Tiny negative eigenvalue может быть
сосчитан и иметь zero residual после budget; исходная матрица сохраняется.

Стоимость light — полный grid scan шести небольших блоков с eigensolves; full дополнительно выделяет
spectral matrix integral, prefix и счётчики, выполняет integral eigensolves. Каждая функция
вызывается явно; observer/checkpoint callbacks и per-event hashing не добавлены. Все возвращаемые
arrays/maps detached: mutation результата не меняет solution/problem/cache.

## Basis subset и preview

```julia
basis_validation_diagnostics(basis::BasisData, projected::AbstractArray{<:Number,3},
    n::NumericalParameters)::BasisValidationDiagnostics
projected_hamiltonian_spectrum(h::AbstractMatrix{<:Number})::ProjectedHamiltonianSpectrum
```

Basis helper принимает `(Nk,Nb,Nb)` projected Hamiltonian и `(Nb,Nb)` T±, возвращает detached
`report::ConvergenceReport` и `limits::Dict{Symbol,Float64}`. Измеряет subset существующего
`validate_problem`: basis orthogonality/eigen residuals, T-adjoint, Hamiltonian Hermiticity и только
для PzP — translation/nearest/two-pairs. Structural limits 1e-12 и norm floors 1e-14 прежние.
Translation limit 0.1; two-pairs limit 0.1 при P_basis<3 и 0.005 при P_basis≥3. Полная static
валидация problem этим не заменяется; scattering problem не строится.

Projected spectrum возвращает `eigenvalues`, complex `eigenvectors`, `hermiticity_residual`,
`eigen_residual`, `norm_floor` и `hermiticity_limit`. Gate effective levels прежний: floor1e-14,
Hermiticity limit1e-12, затем `eigen(Hermitian(h))`. Eigen residual равен
`norm(h*V-V*Diagonal(E))/(norm(h)+floor)`. Nonsquare/empty input — `DimensionMismatch`;
nonfinite/non-Hermitian — `DomainError` до eigen, без repair. Eigenvalues используют input energy
coordinates; при preview это scaled E0. Это diagnostic one-particle levels, не interacting peaks.

## Optimized operators и независимый expected

```julia
production_retarded_green(energy::AbstractVector{<:Real}, h::AbstractArray{<:Number,3},
    sigma_R::AbstractArray{<:Number,4}; η::Real=0.0,
    options::ProductionOptions)::ProductionRetardedGreenEvaluation
production_mix_selfenergy!(target::SelfEnergyFamily, candidate::SelfEnergyFamily,
    alpha::Real; options::ProductionOptions)::SelfEnergyFamily
production_lo_selfenergies(problem::NEGFProblem, green::GreenState,
    cache::ProductionCache; options::ProductionOptions)::LOCollisionEvaluation
sampled_linewidth_diagnostics(A::AbstractArray{<:Number,4}, total::SelfEnergyFamily,
    energy::AbstractVector{<:Real}, energy_weights::AbstractVector{<:Real},
    momentum_weights::AbstractVector{<:Real};
    policy::SCBAPhysicsMarkerPolicy)::SampledLinewidthDiagnostics
fermi_occupation(energy::Real, chemical_potential::Real, kBT::Real)::Float64
```

Resolvent вызывает custom production inverse. `G_R` имеет оси `(E,k,a,b)` и inverse input energy
units; `condition_number` и `dyson_scale` — `(E,k)`, scale в input energy coordinates. Это три
свежих arrays. Reference `retarded_green` не подставляется вместо optimized implementation.

Mix меняет только target по прежней linear формуле `(1-alpha)*target+alpha*candidate` и возвращает
тот же target; alpha∈(0,1]. Все три component shapes согласованы. Writable target components попарно
disjoint: R/L, R/G и L/G не должны иметь пересекающуюся storage по `Base.mightalias`. Каждая из
девяти target/candidate component pairs также disjoint. При нарушении этих условий `ArgumentError`
возникает до любой mutation; implicit copies полных tensors не создаются. Candidate read-only,
поэтому внутренний alias его собственных components допускается. Никакой SCBA stop decision,
Anderson/reset evidence или exact-state receipt этим не создаётся. Existing public
`production_selfenergy_residual` даёт raw max-component Frobenius residual с candidate norm +1e-14;
independent scalar literal остаётся на expected стороне fixtures.

LO port использует фактический cache op/q и owning cache source/grid/kernel/shift/algorithm
validation. LO должен быть enabled, Green shape соответствует problem. Сохраняется thermal и
nonthermal dispatch; результат — fresh lesser/greater `(E,k,a,b)` arrays в E0 units. Green/problem
read-only, cache scratch может меняться. Cache не строится, retarded part и charge normalization не
вычисляются. Fixture сравнивает optimized BLAS/literal dispatch с отдельной полной индексной суммой;
нормализованный kernel и supplied shifts являются входами этого oracle, а не independently verified
microscopic physics.

Linewidth использует unchanged sampled-modal owner и возвращает q10/q50/q90, underresolved, blocks,
status и `sampling_method=:spectral_midpoint_quantiles_A_eigenmodes_v1`. Сохраняет все
status/reasons, включая available, invalid_energy_grid, nonfinite_spectrum,
negative_spectral_weight, empty_spectrum, indefinite_spectrum, nonfinite_broadening и
noncausal_broadening. Unavailable значения остаются NaN. Γ/ΔE dimensionless; нет hbar/2π. Marker не
доказывает grid adequacy.

Fermi wrapper требует finite inputs и kBT>0, делегирует прежнему `_fermi` с saturation ±40. Positive
thermal domain соответствует retained LO fixture (70 K) и physical parameter validation. Новый
wrapper не меняет private IEEE behavior при kBT=0, внутренних callers или отдельный
`HotLOPhononBath(T=0)`/Bose population. Новый запрет T=0 для физической модели не вводится.

Фокусированные fixtures находятся в `test/domain/public_options_contract.jl`,
`test/numerics/public_diagnostic_observations.jl`,
`test/numerics/public_optimized_operator_contract.jl`, `test/physics/public_basis_diagnostics.jl`.
Expected стороны: scalar inverse/quadratic causal root, literal 2×2 adjugate, noncommuting
Keldysh/eigenvalue formulas, full-index LO sum, analytic Fermi values, Rayleigh bounds и BDD sine
spectrum. Их согласие проверяет operator/observation contract на synthetic inputs; scientific
accuracy, fixed-point convergence, domain adequacy и экспериментальную валидацию не устанавливает.
