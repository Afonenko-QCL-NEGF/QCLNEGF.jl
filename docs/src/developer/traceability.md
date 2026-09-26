```@meta
CurrentModule = QCLNEGF
```

# [Трассируемость «формула → функция»](@id traceability)

Таблица связывает математические определения с реализующими их операторами.
[Архитектура](@ref developer-source-map) описывает владение исходниками.
Наличие связи с формулой не заменяет независимую численную верификацию.
Физический смысл терминов следует читать в [словаре](@ref physics-glossary),
а последовательность преобразований и соответствующие рисунки — в
[главе об операторах](@ref physics-operators) и [визуальном атласе](@ref
physical-atlas).

## Базис и функции Грина

| ID | Формула | Реализация |
|---|---|---|
| EQ-BDD-001 | [гамильтониан БенДэниела—Дьюка](@ref eq-bdd-matrix) | [`build_bdd_hamiltonian`](@ref) |
| EQ-BASIS-001 | [локализация по проекции координаты и ортогонализация Лёвдина](@ref eq-pzp-lowdin) | [`build_localized_basis`](@ref), [`build_basis`](@ref) |
| EQ-HRED-001 | [проецированный гамильтониан](@ref eq-reduced-hamiltonian) | [`project_hamiltonians`](@ref) |
| EQ-SHIFT-001 | [энергетический сдвиг](@ref eq-energy-shift) | [`build_shift_matrix`](@ref) |
| EQ-EMBED-001 | [вложение](@ref eq-embedding) | [`embedding_self_energy`](@ref) |
| EQ-DYSON-001 | [Дайсон](@ref eq-dyson-retarded) | [`retarded_green`](@ref) |
| EQ-KELDYSH-001 | [Келдыш](@ref eq-keldysh) | [`keldysh_green`](@ref) |
| EQ-NUMBER-001 | [фиксированное поверхностное число](@ref eq-number-normalization) | [`number_functional`](@ref), [`normalize_lesser`](@ref) |
| EQ-DENSITY-001 | [координатная плотность](@ref eq-realspace-density) | [`electron_density`](@ref) |

## Рассеяние и собственная энергия

| ID | Формула | Реализация |
|---|---|---|
| EQ-KSCALE-001 | [перемасштабирование ядра](@ref eq-kernel-rescale) | [`build_kernels`](@ref) |
| EQ-KIMP-001 | [примеси](@ref eq-impurity-kernel) | [`impurity_kernel`](@ref) |
| EQ-KIFR-001 | [шероховатость интерфейсов](@ref eq-ifr-kernel) | [`interface_roughness_kernel`](@ref) |
| EQ-KAC-001 | [акустическое рассеяние](@ref eq-acoustic-kernel) | [`acoustic_kernel`](@ref) |
| EQ-KALLOY-001 | [сплав](@ref eq-alloy-kernel) | [`alloy_kernel`](@ref) |
| EQ-KLO-001 | [продольные оптические фононы](@ref eq-lo-kernel) | [`lo_phonon_kernel`](@ref) |
| EQ-CONTRACT-001 | [статическая свёртка](@ref eq-static-contraction) | `QCLNEGF._static_contraction` |
| EQ-CONTRACT-LO-001 | [свёртки LO](@ref eq-lo-contraction) | `QCLNEGF._lo_contraction` |
| EQ-HILBERT-001 | [Гильберт](@ref eq-direct-hilbert) | [`direct_hilbert_transform`](@ref) |
| EQ-SELFENERGY-001 | [тождество запаздывающей компоненты](@ref eq-retarded-selfenergy) | [`retarded_self_energy`](@ref) |

## Фиксированные точки

| ID | Формула | Реализация |
|---|---|---|
| EQ-SCBA-001 | [отображение Якоби SCBA](@ref eq-scba-map) | [`solve_scba`](@ref) |
| CONVERGENCE-001 | [полная политика приёмки и застоя](@ref convergence-policy) | [`scba_convergence_assessment`](@ref), [`poisson_convergence_assessment`](@ref), [`scba_convergence_stagnated`](@ref) |
| EQ-POISSON-001 | [периодическая матрица](@ref eq-poisson-matrix) | [`build_poisson_matrix`](@ref) |
| EQ-POISSON-002 | [окаймлённое решение](@ref eq-bordered-poisson) | [`solve_periodic_poisson`](@ref) |
| EQ-OUTER-001 | [внешнее отображение](@ref eq-outer-map) | [`solve`](@ref) |

## Точные операторы полноразмерного расчёта и параллельное расписание

| Контракт | Формула/семантика | Рабочая реализация |
|---|---|---|
| PROD-SHIFT-001 | [нециклический энергетический сдвиг](@ref eq-energy-shift) | [`build_shift_plan`](@ref), [`apply_energy_shift!`](@ref) |
| PROD-CONTRACT-001 | [полная статическая свёртка](@ref eq-static-contraction) | [`production_static_contraction`](@ref), энергетические блоки `:threads`/`:blas` |
| PROD-HILBERT-001 | [прямая сумма главного значения](@ref eq-direct-hilbert) | [`ProductionFFTHilbertPlan`](@ref), [`production_fft_hilbert_transform`](@ref) |
| PROD-SEED-001 | [фиксированное поверхностное число](@ref eq-number-normalization) | `_seed_green_production`, точные скалярные веса ``q_e`` |
| `PROD-RESIDUAL-001` | [приёмочные невязки](@ref eq-validation-residuals) | [`ProductionResidualWorkspace`](@ref), [`production_residual_suite`](@ref), [`production_selfenergy_residual`](@ref) |
| PROD-SCBA-001 | [отображение Якоби SCBA](@ref eq-scba-map) | [`solve_scba_production`](@ref), кандидаты до смешивания на месте |

Параллельные тесты проверяют не новую физическую модель, а эквивалентное
вычисление тех же операторов на конечной сетке. Редукции невязок с
фиксированным порядком отдельно исключают зависимость приёмки от планировщика.
Полная карта владения буферами исполнителей приведена в [исполнении
полноразмерного расчёта](@ref
production-parallel-execution).

## Наблюдаемые и приёмка

| ID | Формула | Реализация |
|---|---|---|
| EQ-POP-001 | [населённости](@ref eq-sheet-density) | [`sheet_density_matrix`](@ref), [`state_populations`](@ref) |
| EQ-CURRENT-001 | [граничный ток](@ref eq-boundary-flux) | [`boundary_flux`](@ref) |
| EQ-CURRENT-002 | [разрешённый/связевой ток](@ref eq-resolved-current) | [`energy_resolved_current`](@ref), [`bond_current`](@ref) |
| EQ-COLLISION-001 | [баланс столкновений](@ref eq-collision-balance) | [`collision_balance`](@ref) |
| EQ-POWER-001 | [баланс мощности](@ref eq-power-balance) | [`power_balance`](@ref) |
| EQ-LEVELS-001 | [диагностические уровни](@ref eq-effective-levels) | [`effective_levels`](@ref) |
| EQ-VALIDATE-001 | [приёмка](@ref eq-validation-residuals) | [`validate_problem`](@ref), [`validate_solution`](@ref) |

Нормативная формула не повторяется в строке документации: она содержит стабильную
ссылку на её `@id`, а API-страница даёт канонический `@docs` ровно один раз.
