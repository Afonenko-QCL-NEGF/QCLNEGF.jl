```@meta
CurrentModule = QCLNEGF
```

# Операторы и диагностика

Оператор преобразует физические данные согласно определённому математическому
правилу: например, гамильтониан определяет стационарные состояния, а уравнение
Пуассона связывает заряд с электростатической энергией электрона. Эта страница
служит справочником функций после знакомства с [физической моделью](@ref
theory-model). Связи конкретных формул, реализаций и проверок собраны в
[таблице трассируемости](@ref traceability).

## Сетки, гамильтониан и базис

Сетки задают узлы и веса интегрирования по координате, энергии и импульсу.
Гамильтониан БенДэниела—Дьюка учитывает изменение эффективной массы между
слоями. Его собственные состояния образуют исходное подпространство,
из которого строят базис, локализованный относительно координаты роста.
Проекция переводит операторы в этот конечный базис; размер базиса поэтому
является самостоятельным параметром проверки сходимости. См. [матрицу
гамильтониана](@ref eq-bdd-matrix) и [локализацию и ортогонализацию](@ref eq-pzp-lowdin).

```@docs
build_grids
build_profiles
build_bdd_hamiltonian
build_localized_basis
build_basis
project_hamiltonians
```

## Периодичность и функции Грина

В периодической структуре электрическое поле смещает энергетический отсчёт
соседних периодов. Функция Грина описывает распространение возбуждения:
запаздывающая компонента определяет спектр доступных состояний, меньшая — их
заполнение. Спектральная функция связывает эти описания с плотностью состояний.
См. [энергетический сдвиг](@ref eq-energy-shift), [уравнение Дайсона](@ref
eq-dyson-retarded) и [уравнение Келдыша](@ref eq-keldysh).

```@docs
build_shift_matrix
apply_energy_shift
build_shift_plan
apply_energy_shift!
embedding_self_energy
retarded_green
keldysh_green
greater_green
spectral_function
number_functional
normalize_lesser
```

## Рассеяние и уравнение Пуассона

Ядро рассеяния задаёт вклад определённого механизма в собственную энергию
электронов. Свёртка ядра с функцией Грина превращает параметры механизма и
текущее электронное состояние в новый вклад рассеяния. Преобразование Гильберта
восстанавливает дисперсионную часть запаздывающей собственной энергии согласно
принятому причинному соотношению. [Глава о ядрах](@ref theory-kernels)
поясняет физический смысл каждого механизма и порядок индексов.

Уравнение Пуассона описывает другой этап самосогласования: пространственное
перераспределение электронов меняет электростатическую энергию, которая затем
входит в гамильтониан. Периодическое решение требует условия нейтральности и
фиксации произвольной аддитивной постоянной потенциала; см.
[окаймлённую систему](@ref eq-bordered-poisson).

```@docs
build_kernels
build_kernels_production
ScatteringKernelBuild
scattering_backend_evidence_class
build_scattering_kernel_set
impurity_kernel
interface_roughness_kernel
acoustic_kernel
alloy_kernel
lo_phonon_kernel
direct_hilbert_transform
retarded_self_energy
production_static_contraction
fft_hilbert_transform
production_fft_hilbert_transform
production_fft_workspace_bytes
retarded_self_energy_fft
production_residual_suite
production_selfenergy_residual
build_poisson_matrix
solve_periodic_poisson
```

## Наблюдаемые

Наблюдаемые величины получают интегрированием функций Грина с соответствующими
весами и матричными элементами. Населённость относится к выбранному базисному
состоянию, а координатная плотность — к пространственному распределению заряда;
это разные представления одной матрицы плотности. Ток, баланс столкновений и
баланс мощности служат также независимыми проверками решения. См.
[плотность](@ref eq-realspace-density), [ток](@ref eq-boundary-flux) и
[баланс мощности](@ref eq-power-balance).

```@docs
electron_density
sheet_density_matrix
state_populations
effective_levels
boundary_flux
energy_resolved_current
bond_current
collision_balance
power_balance
spectral_maps
```
