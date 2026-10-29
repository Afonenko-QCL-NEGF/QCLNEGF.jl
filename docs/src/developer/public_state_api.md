# P02: in-memory state и observations

Эта страница описывает реализацию [семантического контракта S00](public_state_contract.md).
Нормативные имена, обязательность и типы полей определены в
[public_state_types.jl](../../../src/numerics/optimized/public_state_types.jl),
[production_contracts.jl](../../../src/numerics/optimized/production_contracts.jl) и
[public_state_completed.jl](../../../src/numerics/optimized/public_state_completed.jl). Consumer
должен использовать опубликованные Julia definitions. Ручная копия field/schema map в этом документе
не является вторым источником истины. Совместимость consumers проверяется отдельно по опубликованным
owning definitions и их integration evidence.

Core принимает и возвращает значения в памяти. Формат envelope, I/O, committed history cutoff и
execution/attempt/branch provenance принадлежат Contracts/Runner. `CoreDatum` — ограниченный value
grammar для Core warnings, policy-selection audit и final evidence; это не transport API. Он
сохраняет IEEE NaN/Inf, а отсутствие значения представляет `CoreScalar(nothing)`. Неподдерживаемый
тип значения и совпадающие после преобразования record keys отвергаются.

## Public entry points

```julia
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
    execution::ExactExecutionIdentity, domain_adaptation::DomainAdaptationPolicy,
    kernel_options::Union{Nothing,ProductionKernelOptions})::RestoredSolverState

solve_production(problem; options, production_options, execution,
    checkpoint_sink=nothing, checkpoint_context_sink=nothing,
    observation_sink=nothing, history_observer=nothing)
solve_production(problem; restart::RestoredSolverState,
    checkpoint_sink=nothing, checkpoint_context_sink=nothing,
    observation_sink=nothing, history_observer=nothing)

solve_adaptive_production(problem; options, production_options, execution,
    domain_adaptation, kernel_options, checkpoint_sink=nothing, checkpoint_context_sink=nothing,
    observation_sink=nothing, history_observer=nothing)
solve_adaptive_production(problem; restart::RestoredSolverState,
    checkpoint_sink=nothing, checkpoint_context_sink=nothing,
    observation_sink=nothing, history_observer=nothing)

restore_completed_solution(restart::RestoredSolverState;
    revalidation::Symbol, independent_oracle::Symbol)::NEGFSolution
solution_scientific_assessment(solution::NEGFSolution)
```

`execution` обязателен для exact capture/continuation. При активной адаптации обязателен explicit
`kernel_options`; `nothing` допустим при отключённой адаптации. В restart invocation inputs и legacy
`initial_*` не переопределяются. Validated handle непрозрачен: его создаёт `restore_exact_state`,
numerical entry повторно проверяет binding перед compute. Arrays problem остаются caller-owned
read-only на время solve и restore; их изменение требует новой identity. Legacy fixed-Hartree
`solve_scba_production` сохраняет прежние kwargs. Непроверенные внутренние continuation controls
через него недоступны.

`checkpoint_sink(state::ExactSolverState)` получает detached owning snapshot. Это новый public
payload; legacy consumer, ожидающий `NEGFSolution`, должен перейти на P02 definitions. Альтернатива
— explicit keyword `checkpoint_context_sink::Union{Nothing,Function}=nothing`:

```julia
checkpoint_context_sink(actual_problem::NEGFProblem, state::ExactSolverState)
```

Два non-null checkpoint callbacks одновременно отвергаются до restore/copy/hash/compute. Arity не
угадывается; old callback адаптируется один раз к внутреннему двухаргументному вызову. Оба варианта
требуют explicit execution identity для exact capture; при активной адаптации нужен также explicit
kernel policy. Callback синхронный и обязательный: return value игнорируется, exceptions
распространяются и прерывают solve, не попадают в dropped-observation counter.

`actual_problem` — borrowed actual instance, read-only вместе со всеми вложенными arrays только на
время callback. Он не копируется, не становится owned snapshot и не должен удерживаться или
изменяться consumer. `state` независимо владеет detached arrays; его можно удержать после callback.
Fixed/domain-none и adaptive fresh/restart paths передают problem того домена, которому принадлежат
state/identity. Cold-rebuild barrier относится ещё к старому domain; accepted checkpoints после
rebuild получают новый problem, а не initial input. Domain-final и global-final boundaries не
смешиваются. У global final `actual_problem === returned_solution.problem`; deferred final writer
удерживает только detached state и использует `returned_solution.problem` после возврата того же
solve. Reconstruct/registry lookup или completed-state revalidation для получения context не нужны.

Fresh invocations также поддерживают существующие in-memory seed kwargs `initial_Uᴴ`,
`initial_scba`, optional `cache::ProductionCache` и `consume_initial=false`. Это warm starts для
nearby points; они не обещают exact continuation или сохранение counters/history/mixer trajectory.
По умолчанию caller сохраняет seed buffers; `consume_initial=true` передаёт exclusive ownership.
Cache проверяется на соответствие current problem и не разделяется между активными solves; adaptive
cold rebuild использует fresh cache. Validated `restart=handle` по-прежнему отвергает все input
overrides и legacy `initial_*`, включая overrides cache/consume policy. Callback kwargs не меняют
validated numerical policy. Private saved/resume controls не являются warm-start API.

Scientific callback в этих двух full-solve entry points имеет signature
`history_observer(snapshot::ScientificHistorySnapshot)`, technical callback —
`observation_sink(event::SolverEvent)`. Старые низкоуровневые iteration callbacks не приобретают
гарантию detached ownership автоматически.

## Identity и обязательность state

[Canonical encoder](../../../src/numerics/optimized/public_state_identity.jl) использует tagged
binary values с network byte order и SHA-256. Он сохраняет Float64/Float32 bits, complex components,
array element type, dimensions и logical column-major entries; record fields идут в порядке Julia
definition, dictionary/set values сортируются по canonical key bytes. Unitful quantity содержит unit
string и исходное значение. Нет tolerance, Julia `hash` или decimal roundtrip. Type identity может
консервативно различать разные представления одного математического объекта.

Golden primitive: `Float64(1)` кодируется bytes `663ff0000000000000`; SHA-256 —
`e2ac13292d1245bac0ce7b682e9b0efef04e53d29a40be9d72eb2b6bc21b3533`.
[Fixture](../../../test/domain/public_state_contract.jl) содержит literal bytes и hash, проверки
порядка keys, shape, signed zero и one-bit mismatch.

`ExactRestartIdentity` разделяет full actual problem, numerical policy и execution fingerprints.
Policy включает solver caps/tolerances, algorithms, marker policy, roundoff measurement policy,
adaptation и полный kernel recipe. Callback closures и technical cadence в identity не входят.
`ExactExecutionIdentity` требует explicit Core source/code identity, Julia/BLAS/FFTW identity,
scalar/index types и width/chunk/reduction context. Unknown facts запрещены. Core проверяет
доступные runtime facts (Julia, scalar/index/reduction, BLAS width), а source/library build facts
передаёт trusted caller; Core не читает Git, файлы или host inventory. Fingerprint не удостоверяет
достоверность присланного provenance и не заменяет independent oracle. Положительный worker cap
должен быть уже разрешён и не превышать actual Julia default thread pool; иначе Core отвергает exact
binding вместо скрытого уменьшения numerical width.

`ExactSolverState` конструктор требует полный payload, выполняет одну owning deep copy и считает
integrity digest. Все вложенные arrays, mixer states/residuals, pending raw candidate, warnings и
history отделены от producer. Consumer может менять свой snapshot; после mutation старый integrity
не проходит validation. Параметры нового keyword constructor создают новый unchecked payload;
наличие digest не является научной валидацией или доказательством подлинности.

[Validation](../../../src/numerics/optimized/public_state_validation.jl) сначала сравнивает exact
binding, затем integrity, shapes, mechanism/family order, mixer/history/counter/decision coherence и
rebuild receipt. Missing/unknown/mismatch даёт `ArgumentError` с field path; неверная shape —
`DimensionMismatch`. До rejection нет cache build, numerical map или callbacks. Все mandatory
continuation arrays конечны. Failed terminal postmortem может сохранять nonfinite arrays и
diagnostics; это evidence ошибки, не разрешение продолжить невалидный accepted solve.

На **всех** boundaries owning inner options выводятся из завершённой outer history до текущего inner
stage. Сохранённые strict/diagnostic streaks, SCBA terminal classification и original final receipt
проверяются против этой history. Fresh selected-map measurement не делает противоречивую
историческую классификацию истинной. Pending `next_previous_populations` сравнивается побитно с тем
же `_production_state_populations` accumulation order на accepted Green; finite/shape сами по себе
недостаточны.

На `:final` и `:domain_final` сохранённые original metrics повторно проходят owning
`_validate_solution_metrics` под неизменными strict options. Если strict acceptance не достигнута,
approximate band допустима только при выполненном outer criterion, accepted SCBA и явно включённой
`adaptive_working` diagnostic policy. Проверяются original `report.passed`, status и convergence
flag. Passed approximate report не означает strict convergence даже при `scba.converged=true`; при
подтверждённых strict inner и outer критериях strict-passing metrics также не разрешают объявить
original result unconverged или approximate. Эта проверка historical band не вычисляет fresh map и
не меняет original metrics. Отдельная selected-map revalidation при completed restore может дать
другую оценку и хранится своим receipt.

Обязательный `ExactSolverState.mixer_position::MixerPosition` создаёт owning
[mixer](../../../src/numerics/optimized/production_mixing.jl). Он фиксирует число применённых mixing
steps, ordered retained mix indices, последнюю history-clearing reset position/reason и последний
coefficient fallback. Payload обязан покрывать ровно эти indices. Residual-increase reset очищает
историю **до** append текущего mix; causality-guard reset — **после** текущего mix, поэтому пустая
история после guard допустима и проверяема. Coefficient-solve/bound fallback сам по себе историю не
очищает. Потерянный receipt или payload не восстанавливается empty default. Exact capture старого
unvalidated `resume_scba` требует перехода на validated handle.

## Границы, история и adaptation

| Boundary          | Следующее действие                                                                                                                           |
| ----------------- | -------------------------------------------------------------------------------------------------------------------------------------------- |
| `:scba_accepted`  | Использовать сохранённый stop decision; при pending decision один раз применить saved raw candidate к сохранённому mixer, затем начать `ν+1` |
| `:outer_accepted` | Начать следующий outer шаг на `next_U` либо финализировать уже решённый state                                                                |
| `:domain_final`   | Продолжить pending measured-domain decision либо выполнить committed cold rebuild                                                            |
| `:final`          | Только readonly reconstruction и fresh revalidation; continuation entry отвергает state                                                      |

`published_U` соответствует сохранённым G/Σ/n и gauge. Предсказанный Poisson `next_U` хранится
отдельно; final reconstruction не заменяет им U, не центрирует gauge и не выполняет polish.
Сохранённые counters, strict/diagnostic streaks, pending current/populations, seed diagnostics,
witness-retention indices и policy evaluations не сбрасываются. Старые accepted rows не пишутся
повторно. SCBA rows прежних outer stages передаются через scientific history; checkpoint содержит
полную текущую inner history и outer history текущего domain. Durable history chain проверяет
Runner.

`ScientificHistorySnapshot` содержит mandatory `domain_revision::Int` (actual zero-based expansion),
`discretization_id::String` и `basis_id::String`. Внутри Runner-owned run/execution/point stable row
и witness-update target — composite этих трёх полей и существующих `stage`, `outer_iteration`,
`inner_iteration`. Fixed-domain solve имеет revision 0; adaptive resume сохраняет actual revision,
включая domains, пропущенные preflight. `sequence` — delivery ordinal текущего вызова; он не
заменяет stable target. `:witness_update` может адресовать строку из предыдущего invocation.
Контекст и два existing representation digests вычисляются один раз на actual solved domain, без
дополнительного Hamiltonian/spectral coverage и без hash/copy полного state на каждой строке.
History-only invocation без `checkpoint_sink` получает тот же context.

`DomainAdaptationState` содержит current problem recipe **и** nullable next problem recipe. Каждый
`ProblemRebuildRecipe` объявляет physical/numerical/scales/scattering/models/algorithms и весь
`ProductionKernelOptions`; текущие operator arrays не определяют будущий construction recipe.
`DomainRevision.evidence` имеет ровно один из typed вариантов `DomainPreflightEvidence` или
`DomainMeasuredEvidence`. Их mandatory поля и units заданы в Julia types. Все исходные measured
metrics/status/seconds сохраняются, включая failed IEEE values; estimate не выдается за measurement.

Следующий energy grid точно сверяется с owning proposal по mandatory scalar receipts.
`DomainPreflightEvidence` хранит actual `energy_min_scaled`, `energy_max_scaled` и оба
`required_*_scaled`. `DomainMeasuredEvidence` дополнительно хранит четыре actual energy-tail scalars
(`tail_low`, `tail_high`, `edge_gamma`, `edge_spectral`) и `scba_quality`. Producer берёт их из уже
вычисленных estimate/tails/SCBA; новой физической оценки или изменения criteria здесь нет. Scaled
bounds совпадают с построением recorded physical grid. Старое состояние без этих полей не
дополняется default и не считается exact.

Replay сначала воспроизводит scaled `ceil((upper-lower)/old_step)+1`, затем source перевод bounds в
eV. При measured energy failure он отдельно выполняет declared `growth_fraction` и physical ceil, в
том же порядке операций и с обоими исходными node caps. На целой IEEE границе scaled и physical
counts могут различаться; произвольного tolerance или допуска «±1 node» нет. Current next recipe и
каждый historical successor должны точно совпасть с этим proposal: одно сохранение ΔE не разрешает
недостаточный growth. Private recipe-only geometry guard использует направленное IEEE enclosure
scaled-to-eV roundtrip; full public state требует exact receipt replay.

Committed `:finished`/`:cold_rebuild` сверяется с actual source condition: tail threshold, coverage,
remaining expansion budget и invalid SCBA. Failed `NaN`/`Inf` tail diagnostics сохраняются с
исходной IEEE semantics comparison, без подмены нулём/pass. Invalid или exhausted domain законно
заканчивает execution без применения формул следующего окна. Pending `:assess_measured_domain` ещё
не имеет committed current measurement/decision. Current measured revision согласуется с current
numerical bounds/count, grid/basis/Hartree IDs, SCBA quality и original status/report/strict flag.
Node cap применяется к новым domains; большой исходный domain сам по себе не становится invalid
из-за меньшего будущего adaptation budget. Остальные numerical fields, typed policy и kernel recipe
остаются exact-bound.

Эти gates проверяют сохранённую construction/decision coherence, а не независимую истинность
заявленных historical estimates/tails. Старые массивы не реконструируются ради physics replay.
Scientific acceptance и independent oracle остаются отдельными gates.

Обязательный `elapsed_domain_seconds::Float64` сохраняет измеренную стоимость текущего domain через
момент checkpoint, в секундах. Resume добавляет elapsed нового этапа к сохранённому значению;
pending domain-final не теряет уже выполненный solve cost и не выполняет solve повторно. Это
coverage существующего domain clock, без времени между checkpoint и crash, не вошедшего в
сохранённое свидетельство, и без отдельной стоимости cold builder. Unknown cost не заменяется нулём
при decode. Новый domain после cold rebuild начинает свой clock; старые measured revisions сохраняют
собственные seconds.

Supported checkpoint actions — `:assess_measured_domain`, `:cold_rebuild`, `:finished`. Preflight
estimate выполняется до solved state и не создаёт отдельной численной checkpoint boundary. В
некоммитнутом текущем domain revisions содержат все предыдущие domains; committed cold/finished
добавляет текущую revision ровно один раз. Validation проверяет непрерывность revisions и
соответствие action/boundary/commit flag. Restore pending domain decision не повторяет finished
SCBA/Poisson. Committed rebuild передаёт exact recipe builder, увеличивает logical expansion один
раз и начинает cold domain. Pure rebuild до следующего checkpoint может повториться после crash;
transfer старых Green arrays на новую сетку отсутствует. См.
[adaptive implementation](../../../src/numerics/optimized/public_state_adaptive.jl).

## Completed result и scientific assessment

`revalidation` и `independent_oracle` — обязательные explicit keywords восстановления completed
solution; скрытого default recomputation нет. Например, caller явно задаёт
`revalidation=:selected_map, independent_oracle=:none`.

`restore_completed_solution` возвращает обычный `NEGFSolution`; evaluator остаётся существующий
`solution_scientific_assessment`. В `observables[:restoration_evidence]` лежит typed
`CompletedStateEvidence`: original final evidence, fresh selected-map revalidation и отдельный
independent-oracle receipt. Original report/classification не подменяются fresh report.

`:selected_map` вычисляет сохранённый выбранный physical map на неизменных U/G/Σ/n. Это проверка
согласованности state с выбранными algorithms и existing thresholds, без SCBA/Poisson iteration.
Ошибка измерения имеет explicit error receipt и failed fresh report. Нечисловой terminal payload
доступен для postmortem без derived numerical measurement; original NaN/Inf сохраняются.

`:if_supported` отдельно вызывает квалифицированную independent assembly: сейчас только production,
no scattering, full-cell resolvent, full self-energy, resolved transverse momentum, Kramers–Kronig и
sparse nodal-linear shifts. Dense/conservative path разделяет reference helper и имеет
`:unsupported/:shared_selected_reference_implementation`; finite-chain closure имеет
`:unsupported/:finite_chain_assembly_unavailable`. Остальные combinations —
`:unsupported/:independent_map_not_qualified`. `:none` даёт `:not_measured/:not_requested`. Ни
unsupported, ни error, ни совпадение двух unconverged states не являются scientific pass.
Discretization/experimental evidence этим API не создаётся.

## Observations: ownership, units и missing values

`SolverMetric` и `SolverEvent` состоят из immutable scalars/strings и tuples; arrays/problem/cache в
них не передаются. Nullable metric имеет explicit measurement status/reason. Existing NaN/Inf не
заменяются нулём и маркируются `:error/:non_finite`. Scientific unavailable sentinel сохраняется в
исходной строке с отдельным reason. `sequence` монотонен внутри invocation и общий для новых
scalar/history callbacks; outer/inner context указан явно. Runner добавляет attempt identity.

`ScientificHistorySnapshot` передаёт deep snapshot полной `SCBAIteration`/`OuterIteration`, включая
outer populations/density, witness matrices/eigenvector и marker channels/collisions. Header
содержит classification/streaks, units и missing reasons. `action=:witness_update` обновляет
retained payload старой строки без новой accepted iteration; уже переданное значение не мутирует.
При fan-out consumers копирует Runner. Scientific/checkpoint callback failure виден вызывающему;
technical collector failure увеличивает `dropped_observations` и итоговый `observation_dropped`,
сохраняя численное выполнение. Старые control callbacks сохраняют прежнюю error/ownership семантику.

| Величина                                       | Unit в snapshot                                                           |
| ---------------------------------------------- | ------------------------------------------------------------------------- |
| Hartree U; density n                           | `U_H/E0`; `L0^3 n`                                                        |
| SCBA/outer J; technical event J                | `A/m^2`; explicit `A/cm^2` conversion `10^-4`                             |
| Populations; sheet charges                     | target sheet fractions; `L0^2 n_sheet`                                    |
| Outer `r_U`; остальные dimensionless residuals | `E0`; `1`                                                                 |
| Marker μ/linewidth/ΔE; marker current          | existing explicit `eV`; `A/m^2`                                           |
| Collision/witness payload                      | existing dimensionless solver scales, с сохранёнными status/applicability |

Missing reasons различают `:no_previous_iteration`, `:non_finite`, `:not_sampled`, `:not_measured`,
`:measurement_failed`, `:not_retained_by_policy`. Они не заменяют metric acceptance registry. При
обоих checkpoint callbacks `nothing` нет exact identity/state copy/hash на technical events.
Передача borrowed problem в context callback не добавляет его copy/hash. Полная snapshot стоимость
возникает только на requested checkpoint и учитывается как nonelastic `:exact_state_snapshot` memory
burst; это консервативная dimension estimate, не измеренная RAM.

## Проверки и пределы свидетельств

| Утверждение                                                                       | Свидетельство                                                                                                                                                                                                   | Ограничение                                                                            | Решение                |
| --------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------- | ---------------------- |
| Saved scalar/Anderson state продолжает ту же finite fixture trajectory            | [linear restart](../../../test/numerics/public_exact_restart.jl), [literal Anderson oracle и restart](../../../test/numerics/public_anderson_restart.jl)                                                        | Явные synthetic branch-control tolerances; не физическая точность/сходимость           | сохранить              |
| Domain terminal и global final различаются, cold receipt не повторяется           | [pending decision](../../../test/numerics/public_adaptive_restart.jl), [rebuild barrier](../../../test/numerics/public_adaptive_cold_barrier.jl)                                                                | Малый engineering fixture; инфраструктурный crash/durable commit — Runner              | сохранить              |
| Checkpoint получает actual problem fixed/adaptive/restart, legacy arity сохранена | [context port](../../../test/numerics/public_checkpoint_context.jl)                                                                                                                                             | Borrowed readonly problem живёт только в callback; durable final join проверяет Runner | сохранить              |
| Ownership и invalid state rejection явны                                          | [types/metadata](../../../test/domain/public_state_contract.jl), [snapshot lifetime](../../../test/domain/public_state_snapshot_lifetime.jl), [missing groups](../../../test/domain/public_state_rejections.jl) | Не доказательство authenticity внешнего payload                                        | сохранить              |
| Completed selected map отделён от independent oracle                              | [readonly/map/postmortem fixtures](../../../test/numerics/public_completed_restore.jl)                                                                                                                          | Independent map capability ограничена перечисленными combinations                      | сохранить              |
| Technical collector/mode не меняет fixture trajectory                             | [observer equivalence](../../../test/numerics/public_observation_equivalence.jl)                                                                                                                                | Performance overhead не измерен; scientific callback failure остаётся видимым          | дополнительно измерить |

Новые fixtures не меняют production model/default tolerances и не являются научным опытом. Existing
independent mathematical oracles остаются отдельным источником свидетельств.

## Финальные physical context и stored observations

`final_observation_snapshot` — отдельный in-memory port для явного финального postprocessing. Он не
входит в exact-state graph и не задаёт файловый формат. Определения восьми return types находятся в
[final_observation_types.jl](../../../src/numerics/optimized/final_observation_types.jl), выполнение
— в [final_observations.jl](../../../src/numerics/optimized/final_observations.jl). Contracts и
Runner владеют declaration, FileRefs, записью и publication.

```julia
final_observation_memory_estimate(solution::NEGFSolution, terminal::ExactSolverState;
    maximum_metadata_bytes::Int)::FinalObservationMemoryEstimate

final_observation_snapshot(solution::NEGFSolution, terminal::ExactSolverState;
    options::SolverOptions, production_options::ProductionOptions,
    execution::ExactExecutionIdentity, domain_adaptation::DomainAdaptationPolicy,
    kernel_options::Union{Nothing,ProductionKernelOptions}, source_mode::Symbol,
    maximum_bytes::Int, maximum_metadata_bytes::Int)::FinalObservationSnapshot
```

Все keywords обязательны. Estimate проверяет shapes и checked sizes; scalar LAPACK workspace query
не вычисляет eigensystem и не выделяет matrix. Snapshot принимает только global `:final`, выполняет
один `validate_exact_state`, затем сравнивает `solution.options`, U/n, полный SCBA (Green
condition/scale, все scattering/embedding R/L/G, history/mixer) и outer history побитно с terminal.
Float64/ComplexF64 сравниваются через исходные bits: NaN payload и ±0 различаются. Нет tolerance,
gauge alignment, fresh map, reconstruction missing state или solver iteration.

`source_mode=:produced` связывает report/status/converged, warnings и выбранные original
observations с `final_evidence`. `:completed_restore` требует existing `CompletedStateEvidence`:
original точно равен terminal evidence, current report равен revalidation report. Scalar quality
bands и history проверяются существующим evaluator без нового map; approximate passed report не
становится strict convergence. Original report не заменяется. Error receipt с unknown Hilbert
measurement сохраняет его отличие от display `false`. Measured independent receipt требует
context/report independent reference path; unsupported/error/not_measured не дают научный pass.

Inputs borrowed readonly на весь вызов. Snapshot, все его arrays, metric maps и nested audit/receipt
metadata можно удерживать и изменять: они не alias входы. Immutable options/scalars могут безопасно
разделять storage. Positional typed constructors проверяют structural shapes/links, но сами по себе
не удостоверяют producer, подлинность receipt или физическую корректность состояния. Checked capture
port не заменяет независимый scientific oracle.

`FinalObservationSnapshot` содержит context, views, metrics, positivity, summary, assessment,
nullable revalidation/independent_oracle и explicit issues, а также исходные identity/integrity.
`assessment` — detached existing `solution_scientific_assessment`, пригодный для прежнего display;
его lossy projection не заменяет raw metrics/IEEE/status. Научные gates, итерационная сходимость,
дискретизация и эксперимент остаются разными свидетельствами.

### Physical context, axes и units

| Поля                                                    | Семантика                                                                                                                                                                                                                                    |
| ------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `dimensions`, revision, IDs                             | `(Nz,NE,Nk,Nb)` actual problem; fixed revision0 или actual expansion count. Adaptive IDs берутся из receipt, уже проверенного existing validator; fixed IDs вычисляются один раз. Для IDs нет дополнительного Hamiltonian/coverage solve     |
| `x_scaled,wX`; `energy_scaled,wE`; `momentum_scaled,wk` | Exact actual vectors, включая nonuniform weights. Measures: dx, dε, κdκ/(2π); physical multipliers L0 m, E0 eV, L0^-2 m^-2. Spin не входит в wk                                                                                              |
| `z_nm,energy_eV,k_per_nm`                               | Actual scaled coordinates переводятся соответственно x·L0·1e9, ε·E0+Eref, κ/L0·1e-9. Exact scaled arrays сохраняются отдельно; обратный roundtrip не используется                                                                            |
| `state_index,trusted_energy`                            | Int64 zero-based state coordinates и detached Bool mask[E]                                                                                                                                                                                   |
| `scales,inputs`                                         | Field names задают SI/eV units: E0, L0, J0, field/Vperiod, T_L/T_LO, Eref, donor/fion, z0, period length; `g_s::Int`. Actual Core facts, без восстановления из digest                                                                        |
| `problem_scattering`                                    | Все пять construction flags LO/acoustic/impurity/IFR/alloy. Это не доказательство первоначальных Runner requested flags                                                                                                                      |
| `enabled_scattering,selfenergy_families`                | Первая tuple — ordered five-family projection; вторая — full actual scattering order, включая отдельный :electron_electron. Zero requested kernel может отсутствовать; EE не шестой request flag. Embedding остаётся отдельным state payload |
| `physical_models,algorithms`                            | Existing public options сохраняются целиком, включая LO population, EE applicability и выбранную embedding closure                                                                                                                           |
| Boundary/basis/gauge                                    | Stationary field periodic embedding, finite orthonormal projected basis; localization, energy-shift discretization, absolute energy-reference convention и declared effective display gauge                                                  |

Context проверяет structural counts/shapes/links, а не equilibrium applicability. Он не вводит
нового запрета на raw failed IEEE, limiting temperatures или invalid sampled weights. Ограничения
existing physical constructors остаются прежними; reader проверяет применимость своего oracle
отдельно.

### Stored views и явная стоимость

| Поля / logical axes                             | Значение                                                                                                                                                                                                                                                                             |
| ----------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| Four potentials[z], eV                          | Structure E_c·E0+Eref, external -F(z-z0), accepted Hartree U·E0 и сумма; не predicted next_U                                                                                                                                                                                         |
| Density[z], m^-3                                | Published n/L0³, без скрытой repair/reprojection                                                                                                                                                                                                                                     |
| Spatial-energy density[E,z], eV^-1 m^-3         | Один existing full-coherence `spatial_energy_density`; energy measure E0·wE. При согласованном published density интеграл даёт electron density                                                                                                                                      |
| Spectral/occupied traces[E,k], eV^-1            | ReTr(A)/E0 и ReTr(-iGL)/E0; без spin, quadrature weights или 2π. Local `spectral_maps(z,E)` не подстановка                                                                                                                                                                           |
| Localized/effective levels[state], eV           | diag(H(k0)) и one-particle eigenlevels, Eref+E0·value. k0 обозначает первый actual momentum sample (обычно κ=0); actual κ сохраняется в context. Один `project_hamiltonians`, один legacy `eigen(Hermitian(h0))`, upper-triangle convention; новый Hermiticity threshold не вводится |
| Complex wavefunctions[z,state], nm^-1/2         | χ/√(L0·1e9), χV/√(L0·1e9). Phase fixing до unit conversion: first argmax spatial magnitude, positive real pivot. Ascending energy, без matching degenerate subspaces между points                                                                                                    |
| Complex sheet matrices[state,state_prime], m^-2 | Full coherence GL integral и V†PV; populations — real diagonal, fraction denominator Ndop²D·fion, не сумма selected states                                                                                                                                                           |
| Current/sheet scalars                           | e·boundary_flux(embedding_plus), A/m², positive electron flow toward +z; real trace sheet matrix, m^-2                                                                                                                                                                               |
| Positivity                                      | Ровно один `positivity_diagnostics(...;integrate=true)`: все six raw witnesses/statuses/limits/counts и четыре full integral arrays. k→E scan и strict `>` first-on-tie сохранены                                                                                                    |

Outputs передаются owning без второго deepcopy. Новых G/Σ/K или retained witness copies нет.
Expected unavailable views имеют `nothing` и конкретный `FinalObservationIssue` с reason/detail.
Полученные NaN/Inf samples сохраняются; error sentinel не становится measured zero/pass. Unexpected
bugs, interrupt, allocation errors остаются exceptions. Completed restore не запускается этим port:
вызывающий сначала явно выбирает разрешённую runtime operation и приносит result/receipt.

### Checked memory forecast и buffer inventory

Пусть z,e,k,b — actual dimensions. Все products/additions выполняются checked Int до dominant
allocation. Estimate возвращает отдельно caller residency, C/V/D/M/W/B и peak additional allowance.
Caller residency — `Base.summarysize((solution,terminal))`, без повторного owning state copy.

| Allowance            | Numeric payload / owning source                                                                                                                                                                                                |
| -------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| C                    | 8(3z+3e+3k+b)+max(e,8ceil(e/64)); nine Float64 coordinate/weight arrays, Int64 states, word-padded BitVector allowance                                                                                                         |
| V                    | 8(5z+ez+2ek+2b)+16(2zb+3b²)+16; potentials/density, E-z density, two traces/level sets, two wavefunction sets, transform/two sheet matrices, two scalar values                                                                 |
| D                    | 16kb²+8(kb+ek+6k); four positivity integral arrays. Six-kind count/witness records входят в M                                                                                                                                  |
| A                    | 8(z+e+k); existing assessment создаёт fresh quadrature copies, передаваемые без повторного copy                                                                                                                                |
| W baseline           | 16(kb²+8zb+64b²)+8(4z+4e+4k+64b)+65536; full projected H, projected Φ/χV intermediates, h0 copies, local PSD/density/current blocks, conversion/phase/ordering vectors                                                         |
| Q, дополнительно в W | Actual queried WORK/RWORK/IWORK plus W/support arrays, maximum для eigen и eigvals. Eigenvector output уже учтён V; input matrix copies — baseline W                                                                           |
| M                    | 65536+8·summarysize(selected source metadata); report/options, selected audit/receipts, original evidence/warnings и model flags. Проверяется `maximum_metadata_bytes`; headers/strings/records не выдаются за numeric payload |
| B                    | 4·summarysize((problem,terminal))+32768; existing canonical identity/integrity buffers. Это реальная baseline стоимость, не O(1) и не новая streaming-hash implementation                                                      |

Поле W = W baseline + Q; peak additional = C+V+D+A+M+W+B. Guard `maximum_bytes` проверяется до
canonical validation, Hamiltonian, projections, dominant matrices и PSD. Библиотечный workspace
query использует только scalar buffers и -1 workspace lengths; dimensions/products проверяются на
representability в actual `BlasInt` до FFI. Julia1.13 Hermitian eigen/eigvals default driver не
меняется. Query-only return до matrix access предусмотрен
[LAPACK ZHEEVR](https://github.com/Reference-LAPACK/lapack/blob/v3.12.0/SRC/zheevr.f).

Inventory по owning source: `project_hamiltonians` создаёт H[k,b,b], projected/temporary operators;
`spatial_energy_density` — output[E,z], integrated[b,b] и returned coordinate vectors;
`sheet_density_matrix` — complex block/integral и physical conversion matrices; integrated
positivity — четыре output arrays и bounded block-local six-kind/eigvals scratch; boundary current —
energy ordering и local block products. W сохраняет conservative baseline даже после добавления
actual Q. Existing representation-digest buffers учитываются в baseline allowances; нет
дополнительного Hamiltonian ради representation IDs и нет full input-array deepcopy.

Численные live buffers не являются total RSS: compiler, allocator/GC, shared BLAS caches и runtime
reserves учитывает вызывающий resource owner. Core не probe host и не вводит campaign limit. Runner
добавляет собственную wire metadata и bounded serialization slab; файловый commit здесь не
выполняется.

[Constructor/binding/budget controls](../../../test/domain/final_observation_contract.jl) и
[literal full-coherence/units/IEEE/receipts/context controls](../../../test/numerics/final_observation_snapshot.jl)
проверяют этот port без SCBA/Poisson и без completed selected-map evaluation. Искусственные
strict/approximate reports проверяют classification только; они не numerical accuracy evidence.
Publication и downstream native interoperability требуют отдельной owning проверки.
