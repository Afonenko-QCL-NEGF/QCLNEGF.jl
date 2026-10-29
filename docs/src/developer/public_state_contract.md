# S00: публичный in-memory state, восстановление и observations

P02 реализует этот sketch в [public typed API](public_state_api.md). Его реальные Julia definitions
задают обязательные поля и поддерживаемые checkpoint variants; consumer READY по-прежнему требует
owning review и публикации реализации.

## Исторический S00 baseline (2026-10-27)

Дальнейший текст сохраняет утверждённый design-to-types sketch и его исходные source/evidence
ограничения на момент S00. Формулировки «не реализованный» и «нового numerical/runtime evidence нет»
относятся к этому историческому этапу. Текущие definitions, реализация и её проверенные ограничения
описаны в [public typed API](public_state_api.md); результаты проверок P02 фиксируются отдельным
отчётом owning реализации.

Статус: **semantic sketch для P02 design-to-types**, не реализованный и не frozen consumer API.
Source baseline: Core `0df7dd5a481c2c339e8520d2291c035253e7c0e2`, Runner
`8264f48c01ccc524548036d4b3f44552a76f19b5`; root `69d841fdd0fa84f835803560faf80bcfd682679c`. Нового
numerical/runtime evidence нет. Этот документ задаёт предложение Core. P02 определяет и публикует
полные public declarations; только затем P01/P04/P05 получают consumer READY. Файловая оболочка
определяется Contracts после публикации P02, без угадывания вложенных полей/конструкторов по этому
sketch.

## Граница и минимальные signatures

Core владеет in-memory state, восстановлением control flow и численной оценкой. Runner владеет
decode/encode, файлами, history locators, execution/point/attempt/branch identity и durable commit.
Core не импортирует Runner, Contracts, YAML/HDF5, scheduler или host probes. Новая точка входа не
экспортирует private helpers. Existing `NEGFProblem`, `SCBAResult`, `NEGFSolution`,
`solution_scientific_assessment` и математический reference-путь сохраняют свои роли.

Предлагаемые signatures для реализации/уточнения P02 (не frozen declarations и не исполняемый
пример):

```text
exact_restart_identity(problem::NEGFProblem, options::SolverOptions,
    production_options::ProductionOptions, execution::ExactExecutionIdentity;
    domain_adaptation::DomainAdaptationPolicy,
    kernel_options::Union{Nothing,ProductionKernelOptions})::ExactRestartIdentity

validate_exact_state(problem::NEGFProblem, state::ExactSolverState;
    options::SolverOptions, production_options::ProductionOptions,
    execution::ExactExecutionIdentity, domain_adaptation::DomainAdaptationPolicy,
    kernel_options::Union{Nothing,ProductionKernelOptions})::Nothing

restore_exact_state(problem::NEGFProblem, state::ExactSolverState;
    options::SolverOptions, production_options::ProductionOptions,
    execution::ExactExecutionIdentity,
    domain_adaptation::DomainAdaptationPolicy,
    kernel_options::Union{Nothing,ProductionKernelOptions})::RestoredSolverState

solve_production(problem::NEGFProblem; restart::RestoredSolverState,
    checkpoint_sink=nothing, observation_sink=nothing, history_observer=nothing)::NEGFSolution

solve_adaptive_production(problem::NEGFProblem; restart::RestoredSolverState,
    checkpoint_sink=nothing, observation_sink=nothing, history_observer=nothing)::NEGFSolution

restore_completed_solution(restart::RestoredSolverState;
    revalidation::Symbol, independent_oracle::Symbol)::NEGFSolution

solution_scientific_assessment(solution::NEGFSolution)::Dict{String,Any}
```

Последняя функция уже публичная. Нового параллельного scientific evaluator нет.
`RestoredSolverState` — opaque Core handle: включает проверенные problem/options/context и
собственную копию state; Runner передаёт handle обратно Core, не читает его private поля. Оба solve
entry с `restart=...` не принимают одновременно старые `initial_*`, `resume_scba` или несовпадающие
options/recipe. Новый формат не поддерживает silent warm-start fallback.
`restore_completed_solution` принимает только global boundary `:final`. Fixed-domain
`solve_production` принимает `:scba_accepted`/`:outer_accepted` при отключённом adaptation;
`solve_adaptive_production` принимает эти же boundaries внутри active domain **и** `:domain_final` с
pending adaptive action. Handle содержит текущий problem и явный rebuild recipe. Проверка
state/создание handle не выполняют SCBA/Poisson, rebuild или final revalidation.

`checkpoint_sink(state::ExactSolverState)` вызывается синхронно только на границах ниже и передаёт
detached snapshot. Обычный `NEGFSolution` промежуточного checkpoint не содержит всех continuation
полей, поэтому общего `capture_exact_state(solution)` не предлагается. Core собирает snapshot там,
где фактически доступны control state, accepted arrays и next action. Внутренний final активного
adaptive wrapper публикуется как `:domain_final`; `:final` означает завершение всей execution.

## Полный state и accepted boundary

`ExactSolverState` — immutable aggregate с **snapshot-owned mutable arrays**. Immutable struct не
делает вложенные Arrays/Dict immutable. Следующая таблица определяет обязательные public поля;
указанные existing domain types сохраняют полный состав, а не выборку полей для удобства формата.
Таблица фиксирует semantic inventory; полные типы новых вложенных records и public keyword
constructors — обязательный выход P02, см. gate ниже.

| Поле                         | Тип / точный смысл                                                                                                                                                                                                |
| ---------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `identity`                   | `ExactRestartIdentity`, см. следующий раздел                                                                                                                                                                      |
| `boundary`                   | `Symbol`: только `:scba_accepted`, `:outer_accepted`, `:domain_final`, `:final`                                                                                                                                   |
| `published_U`, `published_n` | `Vector{Float64}`, длина `N_z`; поле Hartree и плотность именно сохранённого Green state                                                                                                                          |
| `scba`                       | Полный `SCBAResult`: `GreenState` со всеми `Gᴿ/Gˡ/Gᵍ/A`, `condition_number`, `dyson_scale`; все scattering families, embedding total/plus/minus, history, classification, restart policy, mixer                   |
| `outer_history`              | `Vector{OuterIteration}` текущего domain, без пропуска `populations` и `density`                                                                                                                                  |
| `progress`                   | `SolverProgress`: `completed_outer`, `current_outer`, `accepted_scba`, `outer_strict_streak`, `inner_strict_streak`, `diagnostic_streak` — `Int`                                                                  |
| `continuation`               | `Union{SCBAContinuation,OuterContinuation,Nothing}`, строго соответствует boundary                                                                                                                                |
| `status`                     | Фактический solver status; сам по себе не выбирает boundary и не является scientific pass                                                                                                                         |
| `warnings`                   | Detached records всех уже выданных Core warnings; не technical trace                                                                                                                                              |
| `adaptation`                 | `Union{Nothing,DomainAdaptationState}`, полный declared adaptation state, см. ниже                                                                                                                                |
| `measurement_context`        | `AcceptanceContext`: executed backend/algorithms, Hilbert roundoff measurement flag, map-sharing и U/G alignment provenance                                                                                       |
| `final_evidence`             | `Union{Nothing,FinalEvidenceSnapshot}`; обязательно для `:domain_final`/`:final`: original `ConvergenceReport`, original classification и `final_fixed_point_audit`, charge/model diagnostics с их status/reasons |

Shapes: Green и каждый R/L/G component **каждой** self-energy family имеют `(N_E,N_k,N_b,N_b)`;
Dyson diagnostics — `(N_E,N_k)`; outer populations — `(N_b,)`, density — `(N_z,)`. Механизмы
совпадают с `problem.kernels.enabled` без лишних/пропущенных имён. Mixer хранит `method`, ordered
family names (`enabled...`, plus, minus), одинаковое число `states/residuals`, все family components
и полный Anderson history до принятой границы. Empty linear mixer history допустим; missing Anderson
history не заменяется linear/empty default.

### `:scba_accepted`

Граница после вычисления raw residuals/history, marker policy и `SCBAStopDecision`, **до** terminal
return или in-place mixing. `published_U` — fixed Hartree текущего SCBA, Green и mixed self-energies
согласованы с этой оценкой. `completed_outer` не включает текущий outer шаг; `accepted_scba` —
последний `ν`, а не возраст процесса. `SCBAContinuation` содержит:

- `decision::SCBAStopDecision` целиком, с проверкой против policy/history;
- `candidate::Union{Nothing,FixedPointAuditCandidate}`: уже вычисленные fresh scattering, embedding
  total/plus/minus и roundoff; обязателен, если `decision.terminal == false`;
- `next_previous_current::Float64`, `next_previous_populations::Vector{Float64}` — J/populations
  принятой строки для сравнения на следующей SCBA итерации;
- `seed_mu_scaled` и `seed_number_ratio` с явной measured/missing семантикой;
- `witness_retention::NTuple{3,Int}` — first/worst/last indices для существующей retention policy.

Resume сначала применяет уже принятое terminal decision. Если оно terminal, mixing и новый SCBA шаг
не выполняются. Иначе применяется **ровно один** следующий mixing с сохранёнными candidate и
Anderson state, после чего начинается `ν+1`. Повторной записи строки `ν` нет. Сохранение pending
candidate — консервативный выбор sketch: legacy restart пересчитывает raw map; собственный roundtrip
не доказывает bitwise эквивалентность этого пересчёта. Удалять candidate из нового exact state можно
лишь отдельным изменением контракта с независимым доказательством deterministic reconstruction;
memory gain здесь не заявляется.

### `:outer_accepted`

Граница после добавления outer row и решения о convergence, перед переходом к следующему outer
шагу/финализации. `OuterContinuation` содержит `next_U::Vector{Float64}` и `next_action::Symbol`
(`:next_outer` либо `:finalize`). Сохраняются оба поля:

```text
published_U = Uaccepted                       # Green был вычислен здесь
next_U      = (1 - α_P) * Uaccepted + α_P * candidate_Poisson
```

Ни `published_U`, ни `next_U` не подменяются друг другом. Следующий outer шаг начинает новый SCBA на
`next_U`, с предыдущим SCBA как warm state, без повторного расходования старого inner counter. При
`:finalize` итоговый U/G остаётся на `published_U`; нового polish шага нет. Предыдущие
n/J/populations для внешних residuals берутся точно из последней outer row. `outer_strict_streak`
сохраняется явно: одних outer rows недостаточно, поскольку source обновляет streak с дополнительным
условием `scba.converged` для каждого outer шага.

### `:domain_final`, `:final` и assessment

На обеих границах `continuation === nothing` и final evidence обязательно. `:domain_final` означает
завершённый численный solve текущего domain при ещё не завершённом adaptive wrapper; его public
route — `solve_adaptive_production(...; restart=handle)`, см. ниже. `:final` означает завершение
fixed-domain execution либо wrapper после committed решения `:finished`. Terminal, failed и
approximate результаты допустимы как сохранённые состояния. Это не разрешение возобновить exhausted
inner solve или считать intermediate domain финалом всей execution.

`restore_completed_solution(...; revalidation=:selected_map, independent_oracle=:none)` строит
`NEGFSolution` внутри Core, сохраняя U/G/Σ/n и историю. Обязательная fresh revalidation вычисляет
**ту же saved selected physical map и дискретизацию** на сохранённом state с bound
problem/algorithms/production context; затем применяет existing final validation/thresholds. Это
явный вычислительный audit, не SCBA/Poisson итерация и не polishing. Единственный допустимый
`revalidation` первой версии — `:selected_map`; missing/unsupported selected map даёт явный error,
без fallback к другой embedding/scattering model.

Independent oracle задаётся отдельно: `:none` даёт `not_measured` с причиной `not_requested`;
`:if_supported` требует отдельную эквивалентную реализацию saved model/discretization из
опубликованной P02 capability table либо receipt `unsupported` с конкретной
причиной/неподдерживаемым model switch. Unsupported не является pass и не меняет scientific
acceptance. Existing `_reference_audit_candidate` использует legacy bulk embedding и не принимает
все `AlgorithmOptions`; это не универсальный oracle для `finite_chain` и иных selected model
switches. Нельзя заменить ими saved selected map. Independent finite-chain assembly требует
собственного supported oracle/fixture; пока их нет, его independent status — `unsupported`.
Переиспользование production raw-map helper в revalidation допустимо, но явно остаётся shared-kernel
check, не независимым oracle.

State validity не гарантирует прохождение fresh revalidation. Original и fresh receipts хранятся
отдельно; original context/roundoff missing flag не переписывается новой проверкой. Fresh receipt
записывает фактически выполненные backend/map-sharing/algorithms/Hilbert comparison, а не копирует
старый context. Новая revalidation не повышает original `not_measured` до original `pass`.

Construction, final numerical validation и `solution_scientific_assessment(result)` — отдельные
операции. Нельзя выставлять `scba.converged=true` из одного stored `status=:converged`. Core
сопоставляет original classification, историю/strict streak и свежий validation report; противоречие
receipt/history — invalid state, провал свежей проверки — непринятый результат с причиной.
Assessment сохраняет `not_measured`, `not_applicable`, `error`, `pass`, `fail` и existing registry.
`stationary_candidate_accepted` не означает discretization/experimental validation. Existing
`scientific_accepted=false` и `discretization_verified="not_measured"` не повышаются при restore.

`FinalEvidenceSnapshot` не обязан сохранять derived `spectral_maps` и analysis exports: это
постобработка, не вход следующего шага. Mandatory acceptance context/audit/witnesses сохраняются.
Restored `NEGFSolution` получает observables через Core; отсутствие diagnostic evidence не
исправляется нулём или вымышленным measured flag.

## Exact identity, validation и adaptation

`ExactRestartIdentity` имеет `contract_version::UInt16` и три независимых immutable fingerprints:
`problem`, `policy`, `execution` (`NTuple{32,UInt8}`, SHA-256 canonical binary values). Нет Julia
`hash`, YAML roundtrip, decimal rounding или tolerance-based identity.

- `problem`: **все** поля/reachable values actual `NEGFProblem`: physical/scattering/models, scales,
  numerical parameters, grids/weights/masks (включая angular grid), material profiles, полная basis
  data, ordered enabled kernels и все operator values, четыре shift maps и их discretization.
  Logical matrix axes, dimensions, element types и IEEE bits обязательны; dense/sparse storage может
  иметь один fingerprint только при одинаковых логических entries.
- `policy`: все `SolverOptions` (включая max counters), `AlgorithmOptions`, marker policy,
  `verify_fft_roundoff`, declared `DomainAdaptationPolicy` и весь `ProductionKernelOptions` recipe
  при возможном rebuild. Effective inner tolerance band выводится из immutable enclosing policy и
  сохранённой outer history, не подменяет её identity.
- `execution`: explicit `ExactExecutionIdentity` — Core source/code identity, Julia и numerical
  library identities, scalar/index types, backend, effective worker/BLAS widths и chunk/reduction
  settings. Runner передаёт source/runtime facts явно; Core не вызывает Git и не читает host/files.
  Unknown execution identity отвергает exact restore. Root/Runner packaging SHA не меняет model
  identity; commit Core и numerical execution context квалифицируют trajectory equivalence.

Fingerprint вычисляется при связывании problem/policy, не при каждом technical event. Arrays problem
остаются caller-owned read-only на время solve; их mutation требует новой identity. Смена
grid/basis/window/operators/policy/caps/context — явный mismatch, а не округлённое сравнение. Новый
budget/иной алгоритм требует явной новой постановки/attempt; лимиты не сбрасываются при restart.
Lossless payload означает сохранение значений state без квантования. Bitwise равенство resumed
trajectory — отдельный claim, требующий deterministic context и fixture; до P02 оно не доказано.

`DomainAdaptationState` содержит `policy::DomainAdaptationPolicy`,
`kernel_options::ProductionKernelOptions`, `expansions::Int`, current-domain identity, полный
ordered `revisions` snapshot и nullable `next_problem_recipe`. Recipe содержит explicit
physical/numerical parameters, scales, scattering/models, algorithms и kernel construction policy
для следующего domain. `ProductionKernelOptions` сохраняется целиком: `initial_lookup_nodes`,
`maximum_lookup_nodes`, `relative_tolerance`, `lookup_power`, `validation_fractions`,
`angular_validation_pairs`, `impurity_angular_tolerance`, `strict`. Existing current operators не
определяют эту recipe: разные construction tolerances могут дать одинаковые current arrays, но
разные следующие operators. Required keyword `kernel_options=nothing` допускается только при
отключённом adaptation и отсутствии будущего rebuild, с явной not-applicable семантикой; default
reconstruction запрещена. Builder следующего domain получает буквально bound recipe, в том числе все
kernel options.

`DomainRevision` сохраняет поля, создаваемые `domain_adaptation.jl`: revision/source,
discretization/basis/Hartree identity, actual bounds/NE, required bounds для estimate,
status/metrics/current/momentum coverage для measured state, reason/seed-transfer и original seconds
при их наличии. Optional поля имеют explicit absence; состав определяется вариантом
`:preflight_estimate`/`:measured_state`, не arbitrary extra keys. Next adaptation action
(`:estimate`, `:assess_measured_domain`, `:cold_rebuild`, `:finished`) и decision identity/commit
position хранятся явно. На `:domain_final` pending `:assess_measured_domain` продолжает wrapper с
coverage/decision, **не повторяя законченный SCBA/Poisson domain**. После committed decision
`:cold_rebuild` сохраняются target recipe и revision; restore не добавляет revision второй раз и не
повторяет committed measured action. После successful rebuild expansion counter увеличивается один
раз и новый domain начинает cold seed. Если builder output ещё не вошёл в committed checkpoint,
чистый build той же recipe может повториться после crash; это не новый scientific domain solve и не
повод второй раз увеличить logical counter. Committed `:finished` публикует global `:final`; тогда
adaptive entry не продолжает execution. Checkpoint tests должны покрыть crash до/после decision и на
rebuild barrier. Existing changed domain остаётся cold rebuild, никогда exact transfer старых arrays
на новую сетку. При активном adaptation missing policy/kernel recipe/revisions/position — failure,
не `expansions=0`.

Validation проверяет version, identity, completeness, shapes/order, finite mandatory numerical
arrays (включая **все** Σ и mixer arrays), границу/counters/decision и размеры history vectors.
Строки SCBA/outer contiguous в пределах owning domain/outer stage; source histories привязаны к
explicit outer counter. SCBA истории прежних outer шагов входят в scientific history; их не нужно
копировать в каждый checkpoint. Core проверяет exact local state; Runner отдельно валидирует
parent/attempt/branch и committed history cutoff без directory glob.

Missing/unknown/identity mismatch/invalid arrays — `ArgumentError` с field path и reason;
несовпадающий размер — `DimensionMismatch` с field path. Не происходит seed, solve, file mutation
или quarantine. Нельзя валидировать typed unchecked domain constructor как готовый scientific state.
NaN/Inf в историческом diagnostic поле обрабатываются по metric policy, не смешиваются с finite
обязательными arrays; `no_previous_iteration`, `not_measured`, `not_applicable` сохраняются явно.

## Небольшие observations и ownership

Существующий scalar API `SolverMetric`/`SolverEvent` переиспользуется и уточняется в P02:
`SolverMetric.value::Union{Nothing,Float64,Int64,Bool,String}`, `unit::String`,
`measurement_status::Symbol` и nullable `reason::Symbol`; `SolverEvent.metrics` становится
`Tuple{Vararg{SolverMetric}}`. Все остальные event поля — immutable scalar/string values. У event
есть монотонный `sequence::Int` внутри invocation и explicit outer/inner context; Runner добавляет
execution/attempt identity. Event не содержит Array/Dict, Green/Σ/problem/cache, closure или
borrowed view. Consumer может сохранять его без копирования solver arrays.

`history_observer(snapshot::ScientificHistorySnapshot)` получает отдельный typed immutable header
(`stage`, `outer_iteration`, `inner_iteration`, `sequence`, current classification/streaks) и
`row::Union{SCBAIteration,OuterIteration}` с **detached** вложенным payload. Это snapshot transfer
одному consumer: он может владеть/менять свои arrays после callback, Core их больше не использует.
Копия делается Core один раз на этой границе, включая outer populations/density, marker
channels/collisions, witness matrices/eigenvector. `copy(row)` или `copy(history)` недостаточны.
Рассылка нескольким consumers — ответственность Runner; private fields Core читать не требуется.

Scientific history не coalesce/drop по technical policy: каждая записанная accepted row и retained
witness policy сохраняются. Compact row может явно ссылаться на retained witness update при смене
first/worst/last; нельзя задним числом менять уже переданный snapshot. Callback failure scientific
history/checkpoint виден и делает persistence outcome explicit; technical collector
failure/off/normal/diagnostic mode не меняет уравнения, шаги, tolerances или acceptance evidence.
`observation_sink(event::SolverEvent)` может быть отключён; технические gaps/drop counters —
visible, не scientific `pass`. Core callback остаётся синхронным; bounded batching/I/O — Runner.

Units сохраняются, не выводятся из имени поля: U — `U_H/E0`; n — `L0^3 n`; Green/Σ и charge —
existing dimensionless scale; populations — доли целевого sheet density; SCBA/outer `J` — `A/m^2`.
UI/technical event `A/cm^2` имеет отдельную явную конверсию `10^-4`. ε — `E0`, radial weights —
`L0^-2` с existing radial measure; spatial weights — `L0`. Witness/marker scales и cadence остаются
existing policy; absence verbose trace не является missing scientific witness. Missing value —
`nothing` + status/reason; измеренный ноль — `0` + measured status.

## P02 design-to-types gate перед consumer READY

S00 фиксирует смысл и ownership; имена/типы в signatures — proposal для owning P02. Этот документ
**не разрешает** P01/P04/P05 угадать точный layout или считать их API frozen. Выход P02 — один
опубликованный Core commit с public exported declarations и keyword constructors, справкой API и
matching fixtures. До него эти consumers остаются не READY; P02 design-to-types готов к работе.

Обязательный опубликованный артефакт P02 задаёт:

- Полные поля/типы/constructors новых `ExactExecutionIdentity`, `ExactRestartIdentity`,
  `ExactSolverState`, `SolverProgress`, обеих Continuation, opaque handle,
  adaptation/recipe/revision, warnings, `AcceptanceContext`, `FinalEvidenceSnapshot`, scalar/history
  observations.
- Все допустимые boundary/status/reason/action variants, conditional fields, explicit `nothing` и
  not-applicable cases; public доступ к **раздельным** original/fresh selected-map/oracle receipts.
- Точные callback signatures и все stage/counter/classification/streak fields scientific header;
  public construction fixture consumer, не читающего private fields и не использующего `Dict{Any}`
  для угадывания nested payload. Изменённый `SolverMetric` получает исчерпывающие status/reason
  rules.
- Версионированный Core-only fingerprint encoder: normative field order, type/variant tags,
  lengths/shapes, fixed byte order, IEEE float/complex bits (включая signed zero), optional/enum
  encoding и negative cases non-finite forbidden inputs. Golden byte/digest vectors входят в Core
  tests. Consumer сохраняет public identity/version; повторять encoder в Contracts не требуется.
- Capability table независимых oracles для точной selected model/discretization, с explicit
  `unsupported` reason для остальных; selected-map revalidation и oracle evidence не смешиваются.

Публикация P02 с этими definitions, relevant bounded fixtures и независимым review — условие
consumer READY. S00 review закрывает semantic design; ни оно, ни writer/reader roundtrip не заменяет
этот gate или scientific acceptance.

## Source anchors и необходимая проверка P02

Primary Core locators: [domain types](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/src/domain/types.jl) (`GreenState`, `SCBAResult`,
`OuterIteration`, `NEGFSolution`), [restart policy](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/src/domain/algorithm_options.jl)
(`_solver_restart_contract`),
[actual boundaries](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/src/numerics/optimized/production_solver.jl) (`solve_scba_production`,
`_solve_production`), [witness retention](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/src/numerics/optimized/production_iteration.jl),
[adaptation](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/src/numerics/optimized/domain_adaptation.jl),
[kernel construction policy](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/src/numerics/optimized/kernels.jl)
(`ProductionKernelOptions`), [selected-map audit](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/src/numerics/optimized/certification.jl)
(`_audit_candidate_for_solution`),
[selected embedding](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/src/numerics/optimized/production_operators.jl)
(`_embedding_self_energy_selected`), [assessment](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/src/numerics/reference/final_audit.jl)
(`solution_scientific_assessment`). Runner baseline consumers:
[loader/metadata](https://github.com/Afonenko-QCL-NEGF/QCLNEGFRunner.jl/blob/8264f48c01ccc524548036d4b3f44552a76f19b5/src/infrastructure/persistence/production_state.jl)
(`_load_production_restart`, `_validate_restart_metadata`) и
[completed construction](https://github.com/Afonenko-QCL-NEGF/QCLNEGFRunner.jl/blob/8264f48c01ccc524548036d4b3f44552a76f19b5/src/composition/production_workflow.jl)
(`_completed_restart_solution`).

P02 должен дать bounded synthetic fixtures без production solve: lossless roundtrip всех state
полей; corruption/missing каждой обязательной группы; one-bit operator/grid/policy mismatch;
linear/Anderson continuation; terminal-before-mix no-extra-step; outer Uaccepted/Unew с ненулевым
различием; outer streak с non-strict предыдущим SCBA; exhausted counters; nested mutation isolation;
disabled/normal/diagnostic equivalence и scientific-history collector failure. Adaptive fixture
останавливается на `:domain_final` до/после measured decision и на cold-rebuild barrier: spy
подтверждает отсутствие повторного domain solve, один committed decision/revision/counter и exact
bound recipe у fake builder. При одинаковых current operators одно bit-изменение любого kernel
policy scalar даёт mismatch; native kernel build для этого не нужен. Нельзя принять только
fixed-domain test. Public-only construction/nested roundtrip и golden fingerprint fixtures закрывают
design-to-types gate; это ещё не scientific run.

Completed-state fixture содержит nonzero-coupling finite-chain case, где legacy bulk и cavity maps
различаются: selected-map spy подтверждает выбранную embedding/scattering map, bound execution
context и отсутствие SCBA/Poisson calls. Independent branch получает `unsupported` до наличия
эквивалентного oracle; supported branch проверяется отдельным literal matrix oracle. Original/fresh
Hilbert measurement flags проверяются раздельно; original missing evidence не повышается до pass.

Независимый oracle не вызывает новый writer/reader или общий production helper как единственное
доказательство: literal scalar recurrence проверяет следующий linear/Anderson step и counters;
manufactured Poisson map проверяет published/next U. Existing
[final alignment oracle](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/test/numerics/final_state_alignment_oracle.jl) и
[literal matrix oracle](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/test/numerics/keldysh_literal_matrix_oracle.jl) остаются
отдельными. [Existing restart fixture](https://github.com/Afonenko-QCL-NEGF/QCLNEGF.jl/blob/8d3b4c8f60df6e4fa08a93f979c1c5e81ad1c2c9/test/numerics/inexact_inner_restart_contract.jl)
использует `≈` и общий solver: это supplementary evidence, не bitwise proof и не научная приёмка.

| Утверждение                                                       | Свидетельство                                                               | Ограничение                                                          | Решение                |
| ----------------------------------------------------------------- | --------------------------------------------------------------------------- | -------------------------------------------------------------------- | ---------------------- |
| Semantic sketch даёт P02 design-to-types вход                     | Actual types/boundaries, explicit adaptive route/recipe и раздельные audits | Correction ожидает повторного independent review; API не опубликован | изменить               |
| P01/P04/P05 могут считать S00 frozen exact API                    | Полных declarations/encoder fixtures ещё нет                                | Consumer READY требует опубликованных definitions P02                | данных недостаточно    |
| U/G alignment, полная history и независимый reference обязательны | Actual accepted-state control flow и existing oracles                       | Runtime новый API не проверен                                        | сохранить              |
| Lossless payload уже доказывает exact resumed trajectory          | Такого evidence нет; pending/control gaps перечислены                       | Нужны P02 deterministic fixtures и independent oracle                | данных недостаточно    |
| Новый checkpoint дешевле по RAM/I/O                               | Измерений нет; pending candidate увеличивает payload                        | Engineering cost надо измерять после реализации                      | дополнительно измерить |
