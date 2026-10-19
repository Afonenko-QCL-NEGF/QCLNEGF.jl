# DEV02 — паспорта независимых oracle O08-01…06

Статус: **CP-DESIGN PASS для bounded scopes O08-01/02/03/04/06**
(2026-10-19, независимый reviewer high). O08-03 принят как complete coupled
functional contract; DEV05 I2 source PASS и GREEN-2 925/925 functional
assertions PASS, independent final CP-CODE/docs CPP pending. O08-05/DEV07 остаются **design_blocked**. Полная готовность шести
oracle и полное завершение DEV02 acceptance **не установлены**.

Уточнение O08-02/04 на Core `a8d2c0de21dfdce3ed3c7c585c555656443f6d89`
заменяет неподдерживаемые carrier dimensions штатными checked builders.
Carrier design SHA256
`17c6dda56a2223bde9d20e98928e659d9e60fb7d6615b72ff94dcf190b829ea1`,
independent review SHA256
`b5b3f58eb879272ed4be2081f36fa5c21c9f3733e961951d4649b6feb2edc348`:
O08-04 accepted; для O08-02 coordinator выбрал требуемую exact-ND
precondition, явно включённую ниже. Это portable owning contract; исходный
review receipt в integration workspace:
`docs/implementation/review-closure-2026-10-19/execution/evidence/qcl-oracle-carrier-design-review.md`.
Text consolidation SHA256
`61508afc44b57954b2d1c5d843c333aa8fb4a14fb9a1b1ffa811afc0798fece3`
получила **PASS TEXT/DESIGN**; уточнение conditioning ratio O08-06 получило
**CP-DESIGN PASS**, independent review SHA256
`374672b3d4cce86fd8f1272c9464d0c72ddc52a89b5f6272300ebdcd06567238`.
Final fixture CP-CODE PASS и factual runtime проверены independent review
SHA256 `2cc83ef336d85cd96d0b9e3dde9b6f677524d8f8dacde8379b4c422f9b459f98`.
Underlying PV functional budget остаётся прежним. Результаты targeted
root execution приведены ниже отдельно от design status.

Reviewed candidate SHA256:
`e3f577a55ca215f4e6c33fdc6883f4f38dfa4413b343c092389aac6e88832cca`.
Независимый CP-DESIGN receipt SHA256:
`80d0281721d77ba425a373824f49225c4713683958c2f646bba99cb1c7f4503e`.
Design status сам по себе не означает исполненный test или scientific
acceptance; текущие bounded runtime receipts перечислены отдельно ниже.
Единственный owning output этого пакета:
`QCLNEGF.jl/docs/src/development/review_oracle_contracts.md`;
navigation — отдельный пакет. Этот документ не меняет production code, tests
или physical contract.

## Срез и границы

Чтение: superproject `e0f244d075dd49874ba4f3847f4efa693f3b3622`,
Core `4af6e94808c6e60affbf1050a5184dbc12083335` (DEV01 docs-only поверх
`3108150…`). Core `git status --short` пустой. Root dirty:
изменены Core gitlink, AiiDA worktree, пять system docs; новые execution
receipts и requirement dispositions. Это параллельная работа координатора,
не доказательство нового immutable root. Submodules прочитаны через Git;
не переключались/не обновлялись. Это исходный historical source receipt.
Carrier locators повторно проверены на `a8d2c0d`; численные операторы не
изменялись. Новая docs consolidation не выполняет imports/tests/solver/Git.

Этап — source/design, L-SOURCE ≤30 min до первого checkpoint, ≤2 design
rounds. Imports, tests, build, solver, benchmarks, environment installation
и subagents не выполнялись. Модель исполнения сохраняется: 0 mV, 70 K,
LO/acoustic/impurity/IFR/alloy, electron–electron=none. Малые algebraic
fixtures ниже — явно manufactured functional inputs, не новый production
model и не steady-state solution. Никакие tolerances не изменяются.

Вход: root AGENTS, architecture; Core CONTRIBUTING; handoff development
DEV02–08, autonomy/P1; qcl-physics-audit и physics-checks; historical
`A=docs/implementation/comprehensive-review/2026-10-14-e48e5b2-01`,
PORTABILITY, RV08 K/S/P/L/J/R09–11/R14/P03, RV09 T01–04,
`fixes/RV09-oracle-gaps.md`. Historical source-only review не runtime pass.
Полный повторный аудит 34 операторов не выполнялся.

## Readiness и dependency graph

| Паспорт | CP-DESIGN status | Готовый bounded scope / блокирующий остаток |
| --- | --- | --- |
| O08-01 | design_accepted | Literal noncommuting Keldysh product, anti-Hermitian/PSD controls |
| O08-02 | design_accepted, exact-ND precondition frozen | Independent sine anchor + constant incompatible source, checked full carrier; assessment routing явно synthetic |
| O08-03 | D2 design_accepted; I2 source PASS, GREEN-2 925/925 functional PASS; independent final CP-CODE/docs CPP pending | Complete coupled literal final U/G/n/candidate/history; healthy global FAIL ровно r_sum; independent negative controls |
| O08-04 | design_accepted | Actual NE3/Nk2 public carrier, signed SI current, literal coherence и genuine Nb1 truncation |
| O08-05 | design_blocked | Четыре absolute channel passports ещё не имеют законченной primary→geometry→SI→tolerance цепочки; LO derivation — локальный кандидат |
| O08-06 | functional/ratio design_accepted, CP-CODE PASS | Literal direct/FFT, center-only continuum/tail и accepted derived ratio interval; targeted GREEN2 passed |

`design_accepted` означает только прошедший независимый CP-DESIGN
bounded functional-regression scope. DEV03/04/06/08 имеют текущую
CP-CODE и targeted runtime evidence в своих envelopes. Runtime acceptance
не следует из design; ниже приведены exact fixture/raw receipt identities.
DEV03 ← O08-01; DEV04 ← O08-02; DEV05 ← O08-03
(D2 CP-DESIGN PASS, I2 source PASS; GREEN-2 925/925 functional PASS, independent final CP-CODE/docs CPP pending); DEV06 ← O08-04; DEV07 ← все четыре O08-05 (blocked, допустим
отдельный accepted channel child packet); DEV08 ← O08-06.

Граф независимости: literal arrays/rational sums → API outputs;
analytic derivative → Poisson source → solve → measured diagnostics →
existing assessment; independent SI literals/coordinate products → public
current; finite-support PV integral → exact omitted-cell/tail error.
`kernel → normalized kernel → same kernel` и `direct helper → FFT` сами
не являются independent RHS. Shared code/formats/builders допускаются как
carriers, но их output не вычисляет ожидаемую правую часть.

### Наблюдение: четыре bounded targeted runs

Root execution, Julia1.13.0, threads1/check-bounds=yes, selected prepared
environment; автор этого status update тесты не повторял. Все четыре
receipts: exit0, raw_capture_complete=true, stop_reason=null. PASS опирается
на reviewed assertion contract и raw assertion output с matching fixture
hash, не только на process exit code.

| Card / owning test file | Current fixture SHA256 | Root attempt / assertions |
| --- | --- | --- |
| DEV03 `test/numerics/keldysh_literal_matrix_oracle.jl` | `6adaeea85abf95e21dd0beb763dada3a73d588cae1ed5a5d536eb0e8841a2a1e` | GREEN2,16/16 |
| DEV04 `test/physics/nonneutral_periodic_poisson_assessment.jl` | `9aab7a4eea3a7fcb93033d2c9f6f3b420f6b1872843168e0deb98461f2602730` | GREEN2,183/183 |
| DEV06 `test/physics/nonzero_public_current_oracle.jl` | `17fed2fae5b45549c65cec576cc108f68922a04e3830c0246313ccafd4ee8d01` | GREEN1,248/248 |
| DEV08 `test/numerics/principal_value_independent_oracle.jl` | `4533b1058eb351d5796df18c86ac3b2e2af767078e16bd4cec25588b8e34f70a` | GREEN2,211/211:13 literal+198 continuum/scalar |

Portable raw evidence identities; receipt.json и original.log соответствуют
указанной card/attempt, original failed attempts сохраняются отдельно:

| Root attempt | Receipt SHA256 | Original log SHA256 |
| --- | --- | --- |
| DEV03 GREEN2 | `2c4a574a116d9300fe64bf9da61436d6eb6d8dd24cb3818a40845b61b29a4ab9` | `5c9f677fea9cc2666d80b36c16392a691c971e8d27728ac6e19e972ed92b938a` |
| DEV04 GREEN2 | `8f69ce737c8fde8797947124797be6fd3c9b65b5a5f2726f322dd48164f986da` | `1c098144d968dd42dd3309324051b5fd10e05fbf841144eb63d031fe8901351b` |
| DEV06 GREEN1 | `b30b51c860de104f11c63a61f824e3cc5c076be0e4e9d05f7616040e7a8b0fa0` | `f6e13cf676efa11178a0ee7bc063d1d6c0425ea642727593fcc76a03c043daf8` |
| DEV08 GREEN2 | `96d71dbe946a481cb5b85d6c13f4bd93ae94030794138d88a809699ba92f2811` | `22639e86f27670899583be70b983683fbf374fff0e2df51e1627615912b0c860` |

Всего658 targeted assertions в четырёх выбранных файлах. Это не full suite,
line/branch coverage, fullCI или scientific acceptance. SCBA/outer scientific
solve и final CI здесь0; S1 и DEV05/07 не получают приёмку из этих runs.

## Общее правило допусков

Здесь сохраняются `SolverTolerances` из `src/domain/types.jl:296` и уже
существующие fixture tolerances. Это project regression/acceptance policy,
не взятая из литературы физическая точность. Нового универсального threshold
не вводится. Неполученная оценка — `not_measured`, не ноль и не pass.

Для короткой независимой суммы допустима явная Float64 оценка:
`u=eps(Float64)/2`, `γ_n=n*u/(1-n*u)`, absolute budget
`γ_n*sum(abs, terms)` с указанным числом elementary rounded operations.
Для O08-01 все произведения/суммы — малые Gaussian integers, точно
представимые в binary64, поэтому ожидается exact entry equality, без rtol.
Для SI/radical/π/FFT cases ниже используются уже существующие absolute +
relative fixture budgets; их нельзя объявлять строгой forward-error bound
для любой FFTW implementation. Near-zero comparison всегда имеет absolute
ветвь, а не деление на почти нулевое ожидаемое значение.

## O08-01 — literal noncommuting Keldysh

**Вопрос и применимость.** Homogeneous stationary relation без отдельного
initial-correlation source: `G< = GR Σ< GR†`, `GA=GR†`; при `-iΣ< ≥0`
получаем `-iG< ≥0`. Оператор сам не проверяет причинность заданного GR.
Owning source: `src/numerics/reference/greens.jl:67::keldysh_green`;
signature `keldysh_green(GR::AbstractArray{<:Number,4},
Σl::AbstractArray{<:Number,4})`. Shape `(1,1,2,2)`, все значения scaled,
dimensionless; E/k/spatial quadratures отсутствуют.

**Независимый RHS.** Вставить буквально

```text
GR = [1-i  1; 0  2-i]
C  = [2    i; -i 1],      Σ< = i C
GR C = [2-3i  2+i; -1-2i  2-i]
B = GR C GR† = [7  3+4i; 3-4i  5]
expected G< = [7i  -4+3i; 4+3i  5i].
```

Это ручное раскрытие восьми scalar products, не вызов target/helper.
`det(C)=1`, `tr(C)=3`, `det(B)=10`, `tr(B)=12`: обе Hermitian matrices
положительно определены по 2×2 principal minors. `GR C != C GR`.
Источник relation: Lee–Wacker `cond-mat/0212059v1`, Eq. (5), и owning
`docs/src/theory/07_greens.md::EQ-KELDYSH-001`; web extraction первичного
v1, PDF page4, явно показывает Eq.(5) и GA=GR†. Downloaded PDF identity
зафиксирован ниже. Числовой RHS следует раскрытию выше независимо от
парсинга текста статьи.

**Assertions/норма.** Entrywise exact equality всех четырёх complex entries;
`G< + G<† == 0`; `-iG<` имеет указанные literal diagonal/determinant и PSD.
Нормы необязательны для exact assertion; supplemental Frobenius discrepancy
без весов должна быть 0. Около нуля exact zero assertion, без floor.

**Controls/mutations.** Healthy literal output; wrong sign `-GRΣ<GR†`
даёт `-expected`, общий factor 2 — `2expected`; missing adjoint
`GRΣ<GR` и reversed order `Σ<GRGR†` обязаны отличаться от retained expected
array. Нельзя менять RHS вместе с mutation. PSD alone не обнаруживает
положительный factor 2, поэтому direct entries обязательны.

**Reuse/граница.** Existing diagonal equilibrium/zero-current anchor
`test/physics/analytic_equilibrium_sum_rule_lo_balance_and_zero_current_anchors.jl`
остаётся control; direct noncommuting assertion — единственный новый gap.
Это verification functional relation, не SCBA/Poisson/experimental validation.

## Checked complete carriers O08-02/04

`NumericalParameters` в `src/domain/types.jl:181–284` требует
`Nz≥3,Nb≥1,P≥1,NE≥3,Nk≥2,Nφ≥4,Nqz≥3`. Нельзя задавать NE2/Nk1 либо
несогласованные metadata ради маленького public current вызова. Все physical,
numerical, scales, grids, profiles, basis, kernels, четыре shift matrices,
discretization и models создаёт `build_problem(validate_static=true)`.
Нельзя вручную восстанавливать `NEGFProblem`, `ModelGrids`, `MaterialProfiles`,
`BasisData`, `KernelSet`, менять builder arrays или отключать static validation.

Три carrier configurations: `(Nz,Nb,localization)=(3,2,:none)` для current,
`(3,1,:none)` для genuine truncation, `(4,4,:real_space)` для Poisson и
synthetic routing. Helper с этими единственными аргументами:

```julia
layers = Layer[Layer(d=(Nz*10.0) * u"nm", material=:manufactured_uniform,
    x_Al=0.0, Eᶜ=0.0u"eV", mᶻᵣ=0.067,
    m_parallelᵣ=0.067, εᵣ=1.0, doped=true)]
p = PhysicalParameters(layers,
    (Nz*1.0e16) * u"m^-2", 1.0, 0.0u"V/m", 0.0u"m", 0.0u"eV",
    70.0u"K", 70.0u"K", 12.9, 10.89, 0.0367u"eV",
    2.0e8u"m^-1", 2.0e8u"m^-1", 1.0e-10u"m", 9.0e-9u"m",
    7.0u"eV", 5317.0u"kg/m^3", 5240.0u"m/s",
    nothing, nothing, 2, typeof(1.0u"m")[])
n = NumericalParameters(N_z=Nz, N_b=Nb, P_basis=1,
    E_min=-0.1u"eV", E_max=0.1u"eV", N_E=3, M_E=0.05u"eV",
    k_max=1.0e8u"m^-1", N_k=2, N_φ=4,
    qz_max=1.0e8u"m^-1", N_qz=3, η_seed=0.001u"eV")
s = ScaleSystem(E₀=0.1u"eV", L₀=10.0u"nm")
scattering = ScatteringOptions(LO=false, acoustic=false,
    impurity=false, IFR=false, alloy=false)
problem = build_problem(physical=p, numerical=n, scattering=scattering,
    scales=s, localization=localization, validate_static=true,
    physical_models=PhysicalModelOptions())
```

Это явно изолированный manufactured functional carrier. Модель S1 сохраняет
0 mV/70 K и все пять LO/acoustic/impurity/IFR/alloy, e-e none; fixture
не меняет её frozen plan, scattering или criteria. Uniform layer εr=1,
whole-cell doping и empty interfaces — заданная manufactured geometry.
Неактивные material constants сохранены из reference parameters;
alloy quantities `nothing` согласованы с disabled alloy, не unknown→zero.

Actual builder outputs проверяются независимо:

| Fields | Retained expectation |
| --- | --- |
| x,wˣ | `(.5:1:Nz-.5)`, `ones(Nz)`, Δx=1 |
| ε,wᴱ,trusted_energy | `[-1,0,1]`, `[1/2,1,1/2]`, `[false,true,false]` |
| κ,wᵏ | `[0,1]`, `[1,3]/(16π)`, annular edges `[0,1/2,1]` |
| qᶻ,wᑫᶻ,φ | `[-1,0,1]`, `[1/2,1,1/2]`, `[0,π/2,π,3π/2]` |
| profiles | Ec0, masses.067, εr1, ND1, xAl0, layer_index1; O08-02 requires exact ND guard below |
| kernels | K/qK/enabled empty; FLO zeros shape `(3,Nb,Nb)` because LO disabled |
| four shifts | bias pair I3; LO nodal pair below; shape3×3 |
| discretization,models | `:nodal_linear`; fixed/parabolic/thermal/e-e none |
| scales | E0=.1eV, L0=1e-8m; J0 and λP independently checked below |

With δ=.0367/.1=.367, LO matrices must match
`W+=[1-δ δ 0;0 1-δ δ;0 0 0]`,
`W-=[0 0 0;δ 1-δ 0;0 δ 1-δ]`. Static validation is a carrier property,
not SCBA convergence or scientific acceptance. Other static fields and
complete Green/result inventories are specified in the applicable passport.

Expected path uses literals, not actual `scales`, `CODATA`, `_C̄`,
`_scaled_physics`, BDD/current/Green helpers:

```text
eSI=1.602176634e-19 C, hbar=1.054571817e-34 J s,
m0=9.1093837139e-31 kg, eps0=8.8541878128e-12 F/m,
E0J=.1*eSI J, E0eV=.1, L0m=1e-8 m, mr=.067,
J0ref=eSI*E0J/(hbar*L0m²), λref=eSI²/(eps0*E0J*L0m),
t=hbar²/(2*m0*mr*E0J*L0m²).
```

Public builder is permitted here only as tiny static no-solve assembly:
maximum Hwin12×12, largest Green array96 complex entries, no microscopic
kernel tensor and no SCBA/outer loop. Existing L-JL counters are shared:
per card ≤2 invocations total, each≤120s/1CPU/2GiB, aggregate≤4min,
16MiB temporary/2MiB logs, one common native slot. No hidden setup, benchmark,
retry or full suite. Preserve failed receipts; source/design is not runtime.

## O08-02 — periodic Poisson: совместимость и assessment

**API.** `solve_periodic_poisson(profiles::MaterialProfiles,
grids::ModelGrids, scales::ScaleSystem, nbar::AbstractVector{<:Real})`
и overload `(problem::NEGFProblem,nbar)`;
`src/numerics/reference/poisson.jl:42`.
`U=UH/E0`, `nbar=L0³n`, `x=z/L0`, `wˣ=dz/L0`; positive constant εr,
uniform periodic N≥3, positive donor sheet density. Return `(U,ζ,rP,rneutral)`.

**Independent continuum control — reuse.** Не копировать
`test/numerics/independent_nonlinear_and_continuum_anchors.jl`:
`U*=.01sin(2πx/L)`, `n=ND+εr(2π/L)²U*/λP`, N=32/64/128, εr=12.9,
ND=1000. RHS — analytic derivative, не `build_poisson_matrix*U*`.
Existing L2 error ratios >3.9 и gauge/ζ/residual/neutrality assertions
сохраняются. Этот refinement доказывает конкретный constant-ε stencil,
не variable-ε интерфейсы и не outer fixed point.

**Fixed nonneutral source — новый gap.** N=4; `x=[.5,1.5,2.5,3.5]`,
`wˣ=ones(4)`, εr=1, ND=ones(4), Eᶜ=0, mass ratios positive constant;
healthy `n=ones(4)`, incompatible `n=2ones(4)`. Использовать полный checked
`(Nz,Nb,localization)=(4,4,:real_space)` carrier выше. **До любого Poisson
вызова и synthetic result construction** требуется exact fatal precondition:

```julia
problem.profiles.Nᴰ == ones(4) ||
    error("O08-02 fixture requires exact ND=ones(4); stop and return to design")
```

При false сохранить failed attempt и остановить DEV04; не менять n, RHS,
нормировку или допуски на лету. `ND≈1` недостаточно: actual
`ND=(1+δ)ones4`, literal n=ones4 дают `rP=2|δ|/(2|δ|+1e-14)`, то есть
roundoff δ может стать заметной relative residual. Exact precondition —
выбранное coordinator условие reviewer correction. Current DEV04 GREEN2
183/183 на fixture9aab7a4… подтверждает прохождение fatal exact-ND guard
до двух direct5×5 Poisson вызовов и synthetic construction; это witness
конкретного carrier, не обещание для иных scaled inputs. Positive uniform
wX и остальные static assertions обязательны.
Это не изменение production grid/density. Тогда
`L=cyclic(2,-1,-1)`, `sP=4`, `v=ones(4)/2`, `b=λP*ones(4)`.
Из bordered equations и mean-zero gauge независимо следует

```text
healthy:      U=0, ζ=0, rP=0, rneutral=0
incompatible: U=0, ζ=λP/2,
              rP=2/(2+1e-14), rneutral=1.
```

Важно: ζ относится к **scaled bordered system** `L/sP` и равен λP/2,
а не несмасштабированной mean charge. Данные не neutralize.
`λP=e²/(ε0 E0_J L0_m)` с E0=.1 eV, L0=10 nm; при source constants
`λP/2 > SolverTolerances().r_ζ`. Все inputs/denominators положительны.

**Нормы и допуски.** `rP=||LU-b||₂/(||b||₂+λP*1e-14)` — без spatial
weights; `rneutral=|Σwˣ(ND-n)|/ΣwˣND`, без floor; gauge `ΣwˣU=0`.
Проверять healthy gates по действующим t.r_P/t.r_neutral/t.r_ζ;
nonneutral actual metrics превышают соответствующие limits. Для n=2
neutrality exact=1; `rP>1/2` следует analytic RHS с огромным запасом и
не заменяет проектный acceptance threshold. ζ sign/magnitude сверять с
λP/2 с inherited absolute ζ diagnostic budget `1e-10` из existing sine
anchor (разность ζ, не bound самой ненулевой ζ). Bordered matrix имеет
eigenvalues ±1,1,1/2,1/2, condition number2; этот budget многократно больше
Float64 solve roundoff данного малого well-conditioned fixture. Это
regression policy, не новая physical tolerance. При необходимости строгой
LAPACK-specific forward bound — вернуть D.

**Assessment seam и честная область.**
`_validate_solution_metrics(solution::NEGFSolution,metrics;
accept_approximate_scba=false)`, `validation.jl:751`, принимает полный
metric dictionary. Direct Poisson diagnostics записать как
`:r_P,:r_neutral,:r_ζ`, `:r_U=0` только потому, что здесь U=published U=0.
Direct compatible control утверждает только прохождение трёх designated
измеренных Poisson gates. Для route нужен отдельный synthetic typed carrier;
seed/SCBA solve не запускается. Φ=χ=I4 получены real-space builder.
H0=t tridiag4(−1,2,−1), T+=−t e4e1', T−=T+', Z=diag(x), M−1=I4/mr.
Centres=x, spreads0, overlap eigenvalues1, κ1; window_eigenvalues — все
eigenvalues H0, Egrid0=2t. Residual zeros этого builder следуют его exact
identity representation, не приписываются PzP eigenbasis измерениям.

Независимая complete assembly с d=1, `R=2π²d`:

```text
every e,m: GR=−iR I4, G<=iR I4, G>=−iR I4, A=2R I4,
hm=t*tridiag4(−1,2,−1)+t*κm² I4,
ΣR=εe I4−hm−iI4/R, Σ<=iI4/R, Σ>=−iI4/R,
D=iI4/R, condition_number=ones(3,2), dyson_scale=fill(1/R,3,2).
```

All four arrays `(3,2,4,4)` и metadata `(3,2)` проверены finite/shape,
Dyson/Keldysh/Pauli algebra до plain `GreenState`/`SelfEnergyFamily`
construction. Σ — engineered externally supplied counterterm, не actual
embedding/kernel output. Плотность из этого же Green состояния:
`gs/(2π)*(ΣwE=2)*(Σwk=1/(4π))*R=d`, поэтому n=ones4 обоснована;
U=zeros4 совпадает с direct healthy U. Plus/minus делят **все три** Σ
компоненты пополам, embedding — их точная сумма, mutable arrays не alias.

| Остальные required result fields | Значение и смысл |
| --- | --- |
| SCBAResult.scattering | empty dict согласован с явно disabled manufactured kernels |
| SCBAResult.history | empty: SCBA iterations не выполнялись |
| converged/status/quality | `(true,:converged,:strictly_converged)` — synthetic branch input, валидная classification tuple, не measured convergence |
| restart_contract,mixer_state | `nothing`, default linear empty state: restart/iteration отсутствуют |
| NEGFSolution.problem/options | complete checked carrier; `SolverOptions()` с неизменёнными tolerances |
| U,n,outer_history | zeros4,ones4 из algebra выше; empty history, outer iterations не измерены |
| observables | только явный synthetic label, например fixture_kind=synthetic_acceptance_routing; без scientific quality metadata |
| initial report | `ConvergenceReport(false,Dict{Symbol,Float64}(),String[])`, не accepted; seam получает metrics отдельно |
| outer converged/status | false, `:manufactured_routing`; не публиковать/не сохранять как scientific solution |

Other metric values zero — **explicit synthetic pass sentinels**, не
измерения и не `not_measured=pass`. Frozen keys входного словаря перечислить
буквально, не генерировать из production `_solution_validation_limits`, иначе
удаление gate удалит его и из expected input:

```text
r_D,r_A,r_K,r_Σ,r_λ,r_roundoff,r_Jchange_scba,r_population_scba,
r_PSD,r_caus,r_caus_sc,r_caus_embedding,
r_selfenergy_identity,r_selfenergy_identity_sc,r_selfenergy_identity_embedding,
r_antihermitian_G,r_antihermitian_Σ,r_antihermitian_Σsc,r_antihermitian_Σembedding,
r_sum,r_J,r_ImJ,r_roundJ,r_power,r_ImP,r_roundP,
r_P,r_U,r_n,r_neutral,r_ζ,r_Jchange,r_population,
r_tail_low,r_tail_high,r_tail_k,r_edge_gamma,r_edge_spectral.
```

Optional channel/candidate/hot-LO/e-e keys не добавляются. Actual direct
healthy/bad diagnostics override только rP/rneutral/rζ; rU=0 объяснён тем
же U. Test labels явно говорят synthetic routing. Healthy routing table
даёт report.passed=true/no messages **только как проверку routing**.
Вставка actual three bad diagnostics даёт ровно named failures rP/rneutral/rζ,
без unrelated failures. Затем каждый из трёх isolated cases оставляет только
один bad diagnostic, остальные два healthy, и требует ровно его failure.
Valid false SCBA classification control может быть
`(false,:max_iterations,:unresolved)`, не invalid tuple.

Нельзя передавать sentinels в scientific assessment/persistence или считать
этот healthy routing report физическим pass. У engineered Green spectral
sum defect `4π−1`, то есть полный physical validation даже здесь не заявлен.
Same-state local algebra, route semantics и scientific acceptance — разные
утверждения. Эта complete assembly не закрывает O08-03/DEV05.

**Mutations.** Forced mean subtraction нарушает expected nonzero ζ/rP;
dropped neutrality gate обнаруживается по отсутствующему named assessment
failure при retained `rneutral=1`. Для каждого gate выполнить routing
контрпример, где остальные two gates поданы synthetic healthy inputs,
и явно не называть это independent measurement всей модели.

## O08-03 — complete final U/G: D2 CP-DESIGN PASS

**Текущий статус (2026-10-19).** Принят complete manufactured general-h
functional contract, D2 round2/2. Approved D2 SHA256:
`eb9e1a8c39ad08bafa480970392d844c2c9fffdef02d847397ef87c52111bfc9`;
independent CP-DESIGN review SHA256:
`53c3e9b4061a6cdb42239d6d42f5b3f7a09b94eb36315b1ec22bfea103c4f55f`.
Текущий I2 `test/numerics/final_state_alignment_oracle.jl`, SHA256
`81fbb258f86c1760ed2ffb7964e189518fd1ec7608921e42080c969c4a55f852`,
получил source CP-SOURCE PASS, review SHA256
`bb5b5260487eab19e3c0db01141852c02ebe68ea8c991f68e29e59110db9a4f4`.
Executor GREEN-2 дал **925/925 targeted functional assertions PASS**;
подробная observation приведена после frozen D2 body. **Independent final
CP-CODE и docs review pending**. Healthy global solution report намеренно
FAIL ровно по r_sum; test PASS подтверждает именно этот ожидаемый отказ,
а не physical/S1/fullCI/native acceptance.

Полный D1 owning document SHA256
`d51f984c68ac6c5ea9e25ecce21d2f66aca6b635de168f68cc4e711ce9005be3`
сохранён отдельно до замены этого раздела. Failed I1 SHA6d0db8b1… и
original GREEN-1 receipt557f60ed…/raw6918099f… сохраняются; failure не
переписывается successful attempt. D1 original reading doc673bc964…,
authored D1 docd51f984c… и текущий D2 document — разные identities.
Production inputs/S1/physical thresholds не менялись.

Ниже D2 technical §§1–10 перенесены побайтно, кроме увеличения heading
depth на1. Его source-only status/counter/next-step фразы относятся к
моменту D2 freeze до I2 и GREEN-2; текущее состояние задаёт observation
после body. Они не дают новых разрешений на execution. Обе L-JL попытки
теперь использованы2/2; rerun не разрешён. Формулы, strict γ32 cap и
failure inventories D2 не адаптировались по результату GREEN-2.

D2 references к inherited D1 SI literals означают unchanged общий
checked-carrier блок выше: eSI/hbar/m0/eps0/E0/L0/mr. Exact public seams
и unchanged record ABI приведены ниже после D2 body, чтобы passport
оставался самостоятельным; это source bindings, не новая модель/критерий.
Raw GREEN-1 не печатает actual h. Его значение в D2 — независимое source
prediction; GREEN-2 summary вместе с frozen I2 подтверждает прохождение
reviewed exact unit-chain guards, но не является отдельным dump h.

### 1. Что установил GREEN-1 и что он не установил

Core source identity из coordinator brief:
`f5b73369a96d1225bd96fb34f5262c8dcf965f86`.
Frozen I test SHA256:
`6d0db8b15300dce5350bf1bb282cef008da494be57bd96b23780ab1b314049a9`.
GREEN-1: 46 passed / 0 failed / 1 error, fatal exact predicate в test:662
`problem.grids.wˣ == [1.0,1.0,1.0]` false. Ни одного coupled audit не достигнуто.
39.449 s; sampled RSS928108544 B; sampled owned TMP0; raw2183 B complete;
resource guard_error/stop_reason=null. Exit1 — не scientific failure и не
доказательство production defect: нарушена исходная D1 unit-spacing premise.

Receipt SHA256
`557f60ed37c4f325d09d552672730f6ad6db414781a0da84d57946c467b67c95`,
original log SHA256
`6918099fc6b8a06bc96b70424d08104915c64ed0eb98900b9c3b3a466ab7252e`.
Оба файла `/tmp/qcl-DEV05-GREEN-1/{receipt.json,original.log}` read-only.
Исторический D1 SHA96c2c50a… и review11900d84… не переписываются.
Один из исходных двух DEV05 L-JL invocations израсходован; остаётся один.

D1 ошибочно приравнял dimensionless spacing к exactly1 после отдельных
Float64 nm→m преобразований длины слоя и масштаба. Approximate static
carrier assertions это допускают; subsequent exact premise не следовала
из них. D2 не заменяет false predicate на isapprox и не меняет input lengths
ради получения1. Выбирается **тот же public builder с теми же physical и
numerical arguments**, но независимые литералы выводятся для source-derived h.

### 2. Независимая exact source цепочка единиц и новый fatal guard

Owning chain, без вызова target helpers на expected side:

1. `Layer` (`src/domain/types.jl:1–67`) приводит thickness к metre-valued
   Float64 через `_length`; повторная `_length` уже metre quantity не меняет
   coefficient. `_q/_length/_metres` — `domain/constants.jl:63–81`.
2. `ScaleSystem` (`domain/scaling.jl:10–23`) отдельно приводит `10.0nm` к m.
3. `period_length` (`domain/presets.jl:236`) суммирует один layer; для одного
   элемента result — тот же metre coefficient. `scaled_value(:length)`
   делит его на `L0_m` (`scaling.jl:35–40`); `build_grids` делит ещё на Nz3
   (`reference/grids.jl:12–25`).
4. Выбранный locked Unitful1.28.0 (Manifest git-tree
   `57e1b2c9de4bd6f40ecb9de4ac1797b81970d008`): `src/display.jl:8` задаёт
   nano=-9; `pkgdefaults.jl:74`, `user.jl:326–355`, `conversion.jl:5–47`
   задают metre reference / factor `1//1000000000`.
   `conversion.jl:96–102` умножает Float64 coefficient на factor.
   Julia1.13 `base/rational.jl:162–165,187–189` продвигает Rational/Float64
   к Float64; conversion — Float64 numerator/denominator division.
   Это прочитанные локальные primary source files, не исполнение Unitful.

Пусть RN — IEEE binary64 round-to-nearest-even. Полностью независимые
expected expressions из исходных literal arguments:

```text
nm64  = RN(1.0/1000000000.0)
ell64 = RN(30.0*nm64)
L064  = RN(10.0*nm64)
Lp64  = RN(ell64/L064)
h64   = RN(Lp64/3.0)

nm64  = 0x1.12e0be826d695p-30
ell64 = 0x1.01b2b29a4692cp-25 = 3.0000000000000004e-8 m
L064  = 0x1.5798ee2308c3ap-27 = 1e-8 m
Lp64  = 3.0000000000000004
h64   = 0x1.0000000000001p+0 = 1.0000000000000002
```

Python stdlib scalar reproduction подтверждает эту source chain; actual h
не был напечатан GREEN-1, поэтому последняя строка — **source-derived
prediction**, а не приписанный raw log measurement. Её отдельно проверит
новый exact predicate. Дополнительный Julia metadata readout не требуется
и этим D не запускается.

Expected h не читается из `problem.grids`, `basis`, `_scaled_physics`,
`build_poisson_matrix` или solver output. В future fixture вычислить цепочку
выше из literal input arguments **до** public builder, затем потребовать
exact equality всех трёх actual wZ с h64 и всех actual x с
`RN((g-.5)*h64)`. Дополнительно exact physical layer coefficient=ell64,
actual `scales.L0_m=L064`, typed dimensions m, fixed Nz3 и uniform εr1.
Не использовать `h=actual_wZ[1]` как expected data source.

False predicate остаётся fatal STOP с сохранением нового raw output;
никакого fallback на approximate equality, принудительного `fill(1)` или
конструирования unchecked grids. Archived D1 guard и failed test сохраняются;
новый exact guard проверяет заново выведенный spacing и **всю unit chain**.
Это замена falsified hypothesis после D2 review, а не removal of a guard.

`h` далее означает **exact real value Float64 h64** — рациональное число,
точно вкладываемое в BigFloat/Ball. Уравнения выводятся для любого h>0,
но branch/support assertions ниже заморожены для этого exact h и source inputs.

### 3. Единицы и complete checked carrier

Public builder arguments D1 полностью сохраняются: single doped30nm layer,
Ec0/masses.067/εr1/xAl0, sheet `(1e16/64)m^-2`, gs2, ionization1, bias0,
T70K, scales .1eV/10nm; Nz3/Nb1/P1/NE3/Nk2/Nphi4/Nqz3; E-.1…+.1eV,
margin.05eV; kmax=qzmax1e8m^-1; seed.001eV; channels off/default physical
models, localization:none, `validate_static=true`. No S1 dose/mesh/model
or production criterion changes. D2 does not choose a different representation
of 30nm simply to force unit spacing.

Use D1 independent SI constants (e/hbar/m0/eps0/E0/L0/mr), never CODATA or
target `_scaled_physics` on RHS. With dimensionless `x=z/L0`, h=Δz/L0:

```text
t0 = hbar²/(2*m0*mr*E0J*L0²)          # transverse coefficient, E0 units
tZ = t0/h²                            # finite-difference z hopping
lambda = e²/(eps0*E0J*L0)
D = 1/64                              # scaled sheet number, L0^-2
phi = [1/2,1/sqrt(2),1/2]              # discrete Euclidean normalized
chi = phi/sqrt(h)                     # continuous dimensionless envelope
v = [1/4,1/2,1/4] = abs2(phi)
x = h*[1/2,3/2,5/2],        wZ=[h,h,h],     Lgrid=3h
NDbar = D/(3h)*[1,1,1]
nbar* = D/h*v
sum_g h*abs2(chi_g)=1
sum_g h*nbar*_g=D=sum_g h*NDbar_g
```

Physical volume densities are `nbar/L064^3` and `NDbar/L064^3` m^-3;
physical sheet target is D/L064² m^-2. `nbar` is not sheet density and must
not keep D1 unit-spacing array `[D/4,D/2,D/4]`. Spatial quadrature h supplies
the missing length. Physical envelope ψ=χ/sqrt(L0), m^-1/2.

The source target conversion is independently `Nsheet64*(L064*L064)`
with Float64 arithmetic; here it predicts exactly D. Check that direct
physical-input chain, never call `scaled_value` as its own oracle. Actual
ND remains checked uniform positive and close to D/(3h) by the declared
donor construction error; do not overwrite it with ideal values. `dot(wZ,ND)`
and the three-term donor sum have their γ-bound, not guessed zeros.

**Не смешивать две округлённые длины.** Source `Lp64=RN(ell64/L064)`
равен3.0000000000000004. Exact real grid length3h64 и его Float64 sum
не обязаны ему равняться: `RN(3*h64)=3.000000000000001` в этой цепочке.
Поэтому не вводить ещё один false guard `3h64==Lp64`. Grids используют h64;
`_scaled_physics.Lp` / final r_n normalizer используют input Lp64. Последний
равен `D_actual/sqrt(Lp64)`, а не предполагаемому exact D/sqrt(3h).
В manufactured constant-density history разность n равна0, но его единицы
и denominator всё равно должны быть описаны правильно. Это source
rounding discrepancy внутри одного заданного physical input, не разрешение
изменить period, weights или target number.

Также `RN(sqrt(h64))` здесь равен1.0: actual source χ может совпасть с Φ
побайтно, хотя independent real expression есть φ/sqrt(h64). Поэтому
weighted χ normalization и density comparisons имеют propagated division/
sqrt/Φ error, а не новый exact `sum(h*abs2(χ))==1` guard. Ideal relations
в таблице остаются физически выведенными; actual Float64 fields проверяются
со своими направленными roundoff enclosures.

All remaining grid fields are D1 values: E[-1,0,1], wE[.5,1,.5], trusted
[false,true,false], k[0,1], wk[1,3]/(16π), qz[-1,0,1], wq[.5,1,.5], phi-angle
[0,π/2,π,3π/2]; bias maps I3 and LO maps δ=.367. Nodal-linear discretization,
empty kernel dict/enabled and FLO=zeros(3,1,1) are unchanged explained limits.

### 4. Basis/BDD/kinetic terms: h powers must stay distinct

From the literal 3-node BDD cell matrix `tZ*[2 -1 0;-1 2 -1;0 -1 2]`:

```text
Phi = phi up to one global unit phase       # not phi/sqrt(h)
chi = Phi/sqrt(h)
H0 = (2−sqrt(2))*tZ
Tplus=Tminus=tau=−tZ/4
c=2*tau²=tZ²/8
Z=3h/2; centre=3h/2; spread=h/sqrt(2)
Mminus1=1/mr
window_eigenvalue=H0; Egrid0=2tZ
h_m(U)=H0 + phi†diag(U)phi + t0*k_m²
```

**Transverse kinetic coefficient stays t0, not tZ.** Only spatial finite
difference hopping carries h^-2. Changing the k term to tZ would alter the
projected equation and would be a new shared oracle error. `basis.jl:14–84`
explicitly sets χ=Φ/sqrt(Δx) and Euclidean Φ orthonormalization; projection
uses Φ, so no extra h factor multiplies φ†diag(U)φ. `hamiltonian.jl:36–48`
supplies the BDD h^-2; `basis.jl:288–318` supplies the separate k term.

Static basis expectations must use h-scaled χ/Z/centres/spreads and tZ-scaled
H0/T/window/Egrid. φ itself, overlap eigenvalue/condition1, Mminus1 and
translation diagnostics retain D1 meanings. All arrays/types/no-alias and
full model/problem inventories remain required.

### 5. Independent Poisson/Hartree relation

Write literally, never obtain matrix/RHS from target helpers:

```text
M=[2 −1 −1;−1 2 −1;−1 −1 2]
L_h=M/h²
b=lambda*(nbar*−NDbar)=lambda*D/(12h)*[-1,2,-1]
U*=lambda*D*h/36*[-1,2,-1]
sum(U*)=0;  h*sum(U*)=0;  L_h U*=b
zeta*=0; neutrality*=0
u*=phi†diag(U*)phi=lambda*D*h/72
```

This gives one consistent density/target/U/h/Uprojection relation. U is
dimensionless Hartree energy UH/E0; physical eV is U*.1. There is no
new Poisson solve or relaxation in this design; this is literal algebra.

Future arrays `n_float` and U_float are nearest-rounded independent literals
from these high-precision expressions. Once rounded, treat **their Float64
values as exact inputs** to comparison enclosures; do not silently substitute
`sum(n_float)=D/h` or ideal U for the measured final inputs.

#### Actual floating stencil and exact tiny bordered-system reference

The source stencil forms Δx² in Float64. Independent input-side expressions:

```text
h2_64=RN(h64*h64)
qL_64=RN(1.0/h2_64) = 0x1.ffffffffffffcp-1
L_api_literal=qL_64*M
sP_literal=4*qL_64
```

For this qL, multiplying by2/4 is exact and each row absolute-sum is
4qL. Verify this scalar arithmetic from the literals before comparing the
tiny solver; source face εr harmonic mean is exactly1. No call to
`build_poisson_matrix` is needed. `L_api_literal` is an independently
derived **floating discrete operator**, distinct from ideal L_h by a known
roundoff coefficient; it is not a changed physical grid/stencil.

Let `la` and uniform `nd` be **checked actual input coefficients** (bounded
against independent lambda and D/(3h)), let n_f be the literal rounded
published n array, `S=sum_j n_f[j]` evaluated exactly in high precision,
and `P=I−ones(3,3)/3` literal. The exact solution of the actual-coefficient
bordered discrete system, before its floating solve, is independently

```text
b_a=la*(n_f−nd*ones3)
U_a_exact=la/(3*qL_64) * P*n_f
zeta_a_exact=la*(S−3nd)/(4*qL_64*sqrt(3))
neutral_a_exact=abs(S−3nd)/(3nd)
```

Uniform donor errors feed ζ, while U depends on the zero-mean part of n_f.
The h factors in neutrality cancel only in exact algebra; source dot/sum
roundoff is still explicitly included. D1 simplification `la*D/36*d` is
not the D2 actual-coefficient U reference. The ideal `lambda*D*h/36*d`
remains the independent physical/discrete literal state, with rounding and
la/nd/qL perturbations kept in the subsequent enclosure.

For actual returned Up, expand the three rows of `qL_64*M*Up−b_a` independently.
With r this exact residual and measured gauge mean(Up), the forward bound is

```text
||Up−U_a_exact||inf <= ||P*r||2/(3*qL_64) + abs(mean(Up))
```

Add the explicit b formation/row-accumulation error divided by3qL. This uses
the literal positive zero-mean eigenvalue3qL, not guessed LAPACK accuracy.
Compare returned rP against independently expanded residual /
`(norm(b_a)+la*1e-14)`, propagating its floating multiplication/norm error.
Retain unchanged existing rP gate1e-11. For ζ do not assume LAPACK forward
accuracy from epsilon: the source uses `v64=RN(1/RN(sqrt(3.0)))`. Freeze this
independent literal and the 4×4 scaled matrix

```text
Q4=[ M/4, v64*ones3; v64*ones3', 0 ]
rhs4=[ b_a/(4qL_64);0 ]
zeta_exact_v64=sum(b_a)/(12qL_64*v64)
```

Its two zero-mean eigenvalues are3/4 and the other two are±sqrt(3)*v64;
the literal v64 interval proves `norm(inv(Q4),2)=4/3`. Independently form
the exact four-row residual `R4=Q4*[Up;zeta]−rhs4` from the returned floats,
with directed arithmetic, to get `abs(zeta−zeta_exact_v64)≤(4/3)*norm(R4)`.
This measured residual includes source b/rhs formation/solve error and
does not use a second solve. It also gives an independent full-vector
forward enclosure. Reference `zeta_exact_v64` differs from the ideal
sqrt(3) reference by the explicitly computable reciprocal-v error; either
retain the exact v64 reference or add that outward difference.

Retain the D1 short-formation comparison as well, now with its h/q units:
`Bz_short=gamma32*(la*D/(qL_64*h)+abs(zeta_exact_v64))`, outward, including
the explicit ideal-v versus v64 difference when needed. The new measured
forward bound must fit this retained roundoff regression budget; it does
**not** enlarge it. Require `Bz_forward≤Bz_short` and use the smaller enclosure
for comparison. The γ32 budget alone is not advertised as a universal
LAPACK forward guarantee; the measured Q4 residual provides that evidence
for this actual tiny solve. Failure of either retained regression or
forward check is STOP, not permission to widen γ or a physical tolerance.

Enclose neutrality against its exact ratio plus γ32 times the absolute
weighted dot/sum terms divided by positive actual donor sheet; near zero
this is an absolute cancellation budget, not a relative error in zero.
Require the resulting ζ/neutrality upper bounds to be below unchanged
rζ/rneutral gates. A wrong output cannot inflate its residual enclosure
and still receive a functional pass. All exact refs are frozen before the
tiny public Poisson call; residual measurement follows that call and precedes
the independent final validation call.
Enclose rU from `max_j abs(U_a_exact[j]−U_published_float[j])` plus this
forward bound and final subtraction rounding. No new fitted Poisson tolerance.

### 6. Six-node Green state, same-G raw map and charge

Use the same exact branch inventory as D1, now with

```text
z_em=E_e−[(2−sqrt(2))*tZ+u*+t0*k_m²]        # healthy
z0_em=E_e−[(2−sqrt(2))*tZ+t0*k_m²]         # old U=0
c=tZ²/8
g=(z−i*sqrt(4c−z²))/(2c)                 if z²<4c
g=2/(z+sign(z)*sqrt(z²−4c))              otherwise
```

For both exact source-derived h states only e2,m1 lies inside the band;
every other node is real and nonsingular. Branch/radicand/denominator/support
guards stay fatal and are frozen before target calls. With independent
`B=16π²D` at that node:

```text
A=−2Im(g); G<=iB; G>=−i(A−B); f=B/A, 0<f<1
Sigma±X=tau²*X; SigmaTotalX=c*X, X in (GR,G<,G>)
```

At other nodes A/G</G>=0 follows from real causal roots; their nonzero real
GR/ΣR enter all Dyson/Σ norms. Arrays shape(3,2,1,1), no aliases; condition1,
positive Dyson scale `abs(z−c*g)=abs(1/g)`. There is no η-work or imported seed.
At active node `c|g|²=1`, hence raw Keldysh and spectral identities; outside
both sides vanish. Actual fresh map uses `_reference_audit_candidate` /
embedding target on **actual side**; compare every field against independent
τ²G/cG. Empty channels/empty-roundoff aggregate0 and legacy bulk-feedback
scope remain explicit. No microscopic model/scattering oracle is manufactured.

Energy/radial/charge normalization **does not acquire an extra h**:

```text
N(G<)=gs/(2π)*wE2*wk1*B
     =2/(2π)*1*1/(16π)*16π²D=D
P_sheet=[D], population=[1]
nbar_g=D*abs2(chi_g)=D/h*phi_g²
sum_g h*nbar_g=D
raw/candidate/normalized charge=D; lambda_raw=lambda_candidate=1
```

h appears only when reconstructing spatial volume density and integrating
it back with wZ. Do not insert h into `number_functional`, `wE*wk`, current,
sheet matrix or spectral sum. Charge/Pauli do not use fitted normalization.
Predicted numbers for h64 differ from D1 only at last arithmetic digits:
healthy A≈37.5340379643918,f≈.0657376939463093; old A≈55.7172917139956,
f≈.0442842971072220. Norms and scalar floors remain unchanged.

### 7. Required complete records and histories

D1 full-field inventory remains binding; D2 changes only their derived
spatial values. `NEGFProblem` uses the checked builder; every actual field
is bound to the above independent input expectations. GreenState has full
GR/G</G>/A/condition/scale; stored embedding/plus/minus follow τ_h² map;
candidate comes from same published G and is independently checked. All
typed dictionaries/classifications/mixer/restart optional fields keep the
same explained manufactured meanings, never arbitrary zeros.

`NEGFSolution.U` uses independently rounded U*, `.n` rounded nbar*=Dφ²/h,
options unchanged; final/public density and number checks are required.
For `electron_density` physical conversion use canonical L064³, with the
short division/conversion/multiplication rounding in its enclosure. Bound
actual χ²D against intended nbar_g and the rounded literal n_f separately;
do not keep an old D1 density expected array. Sheet conversion uses L064².

Manufactured inner ν1…4 and outer μ1…3 sequences are still constant algebraic
states. First-predecessor differences Inf/NaN; later density/current/population
differences are derived from identical arrays/samples. History chargeD is
sheet number, so stays1/64. Outer density vectors now n_f(h); old-Hartree
outer rU uses `max(abs(U*(h)))`. Retain same-field measured SI J in history,
with independent zero-current cancellation bound (τ now τ_h). This avoids
false near-zero temporal relative change; no history becomes fake solver
runtime. Optional witness/markers absent and provisional report unassessed.

No supplied fake pass metric dictionary. For all variants call actual
`fresh_fixed_point_metrics`, `validate_solution`, and explicit
`_validate_solution_metrics(...;accept_approximate_scba=false)`. The frozen
38-base /41-with-candidate gate inventory and source-limit key comparisons
from accepted inventory delta stay unchanged. Full API outputs supply all
other gates, including known inadequate spectral quadrature.

### 8. Roundoff contract: exact input spacing, same gates

No `SolverTolerances` or inherited ATOL/RTOL changes. New exact unit-chain
guards do not use ATOL/RTOL. D1 static tolerances stay only where already
accepted for eigensystem/constant conversion comparisons, now around the
correct h-dependent reference.

1. Compute h64 using the stated source-equivalent literal Float64 operations;
   verify canonical input coefficients and actual grid **exactly**. Embed its
   exact binary value into 384-bit directed intervals. Do not use decimal
   parse1 as a surrogate, interval centre from actual wZ or a measured h fit.
2. All h-dependent ideal arithmetic is outward: h², 1/h, sqrt(h), tZ=t0/h²,
   U/u=lambda*D*h factors, n=Dφ²/h. Once U/n/g are rounded to Float64, treat
   their actual literal floats as exact inputs. Independent polynomial defect
   of g retains actual construction error; active sqrt cancellation factor
   `(4c+z²)/abs(4c−z²)` stays explicit. No large BLAS/FFT bound is inferred.
3. Constants/root derivation uses directed 384-bit operations, including the
   extra h multiplications/divisions; counting them is mandatory in ledger.
   Float64 target scalar product ceilings remain γ16 per complex triple,
   γ32 per six-entry norm/charge/projection, γ64 per scalar Dyson metric.
   h64 itself is exact data after its exact source-chain guard, not an
   unmeasured perturbation to be absorbed by increasing γ or physical gates.
4. ηh expected Hamiltonian is now H0_h+φ†Uφ+**t0**k². Include actual checked H0,
   Φ error and actual kinetic coefficient/k deviation as D1; expected H0 uses
   tZ, not t0. Input-side actual T coefficients are compared with τ_h and give
   ηΣ. Triple-product, A/K/Σ/charge quotient enclosures retain D1 formulas
   with h-correct g/c/τ. Denominators stay positive, all divisions guarded.
5. Donor construction expectation is D/(3h), with source target unit conversion
   and dot normalization. At most: source L0² multiply and sheet multiply,
   three positive h terms and division — ≤8 elementary rounded operations;
   existing γ16 conversion envelope and γ32 weighted donor-sum envelope
   suffice. Actual sum includes h explicitly. Uniform/finite/positive ND
   guards are required. No `ND≈D/3` unit-spacing shortcut remains.
6. Density/χ expected side is φ/√h and Dφ²/h. Measure actual χ field deviation
   from this independent expression and propagate six-entry integration plus
   physical-unit conversion (γ32). This is not changing intended n to target
   electron-density helper output. Sheet/number measure has no h.
7. Poisson uses independent floating qL_64 M, exact rounded n_f as input,
   actual checked la/nd as perturbation coefficients, three-row residual,
   zero-mean eigenvalue3qL and gauge; formulas in §5 replace all unit-spacing
   exact simplifications. Actual qL formation is frozen independently of
   `build_poisson_matrix`. Residual/norm/source error counts retain γ3 b
   formation and γ32 row/dot/norm envelope, with corresponding 1/(3qL) factor.
   The independently expanded Q4 residual additionally bounds ζ by4/3,
   using the exact source-derived v64; no LAPACK accuracy coefficient is
   guessed and no new physical gate is introduced. Its forward enclosure
   must fit the retained h-correct γ32 short-formation comparison, never
   replace it with a looser measured-error allowance.
8. Current cancellation per side remains
   `2γ3*τ_h²*B*(A−B)/(16π²)` (plus checked T discrepancy as in D1).
   rJ bound divides scaled flux by1e-12, temporal history J floor is
   1e-12 A/m²; they are different units. Bounds remain <1e-6. No exact J0
   or zero-history shortcut is introduced.

Frozen comparison: for any scalar Dyson reference D0,g, n0=abs(D0g−1),
d0=abs(D0)*abs(g)+1, propagate ηh+ηΣ and γ64 to bn/bd exactly as D1;
require d0>bd, bound ratio by
`[max(0,n0−bn)/(d0+bd),(n0+bn)/(d0−bd)]`. Max six intervals. Healthy upper
must be below source gate before pass assertion; wrong/stale lower must
exceed it. A/K/Σ retain full-array Frobenius+1e-14 and their same quotient
logic. Missing/out-of-bound input remains STOP, never tolerance tuning.

### 9. Healthy and wrong inputs, unchanged exact failure inventory

Healthy: U*(h), n*(h), six g_h/Σ_h and actual same-G candidate. All designated
raw/alignment gates pass algebraically; global spectral sum deliberately
fails unchanged1e-3. With only active spectral node,
`r_sum=max(abs(A_h/(2π)−1),1)`≈4.9737276762. No scientific solution pass.

N1 wrong/predicted U: publish2U*(h), retain G_h, Σ_h, n*(h), candidate and
old histories. Fresh Σ residual remains near0; maps depend on same G at
bias0. δh=u*(h), literal rD per node
`abs(u*g)/(abs((1/g−u)*g)+1)`; rU=max|U*|=lambda*D*h/18.

N2 stale G at fixed healthy final U*: construct G0_h/Σ0_h for Uold0 using
H0_h+t0k² and the same B; same nbar*=Dφ²/h, target/charge/population. Publish
healthy U*(h) with this complete old state. rU passes but rD has the same
formula with g0. Stored and fresh rΣ stay small. Old outer records explicitly
retain the predicted unsolved U update; no synthetic accepted history.

For both variants at active node |g|=1/sqrt(c_h), put
`x_h=u*(h)/sqrt(c_h)=lambda*D*sqrt(8)*h³/(72*t0)`.
Triangle inequality gives `rD>=x_h/(2+x_h)`≈.08897096039.
Numerical independent scalar expectations: N1 rD≈.08946567246;
N2≈.09013587908; rU_wrong≈.01570757654490 E0. These are mathematical
separation margins against existing gates, not new physical thresholds.

| Full input | Exact expected global failed-key set |
| --- | --- |
| Healthy U_h*,G_h* | `{r_sum}` |
| N1:2U_h*,G_h* | `{r_D,r_D_candidate,r_U,r_sum}` |
| N2:U_h*,G0_h | `{r_D,r_D_candidate,r_sum}` |

Candidate D limit remains t.rK1e-8; stored D t.rD1e-11; rU t.rU1e-7.
Require fresh and last stored rΣ<unchanged1e-8; exact failed-key sets keep
r_sum from masking missing coupled gates. No gate drops/zero dictionary.
D1 optional N3 remains unimplemented; no extra control or runtime scope is
added in D2. Existing narrow-pole primitive remains separate evidence.

### 10. Exact implementation handoff and stop

D2 supersedes only D1 exact-unit-spacing assumption and all consequences
enumerated in the accompanying narrow math diff. Other D1 field inventories,
actual final API calls, source statuses/norms, model/scattering scope, no-alias
rules, original failure retention and finite resource envelope stay binding.

After independent D2 CP-DESIGN, a medium owner may change only the new
DEV05 test and associated author receipts to materialize this correction;
owning docs changes are a separate ultra authoring step. No source operator
or S1 input changes, no setup/full suite. Preserve I6d0db8…/D1/docd51… and
GREEN-1 as failed originals. At most **one remaining** DEV05 L-JL invocation
≤120s,1CPU/2GiB,TMP16MiB/log2MiB, aggregate original two invocations≤4min;
root dispatch/review owns it. This D uses0. No metadata diagnostic consumes it.

If exact unit-chain guard, branch/support, derived enclosure or final failed
sets do not hold, preserve that final attempt and report the concrete blocker;
no third D round, no automatic rerun or tolerance/grid/model adjustment.
Native/runtime outcome is not inferred from source derivation.

| Утверждение | Свидетельство | Ограничение | Решение |
| --- | --- | --- | --- |
| D1 exact-unit-spacing premise неверна | Actual GREEN-1 fatal test662 and frozen raw hashes | Actual h value not printed; audits never reached | изменить |
| Source-derived h и complete h-dependent contract выведены | Unitful/Julia/core source chain; independent bit arithmetic; χ/n/BDD/Poisson/TGT equations | D2 CP-DESIGN pending; Julia not run | дополнительно измерить |
| Failed originals и owning files сохранены | I6d0db8…, D196c2…, owning docd51… hashes | Publication halted by coordinator | сохранить |
| Healthy fixture/S1 научно принят | Spectral gate still fails; no S1 run | Discretization/experiment absent | данных недостаточно |

### Executor observation — DEV05 GREEN-2, 2026-10-19

Чтение сохранённых executor artifacts, без повторного запуска: matching
fixture81fbb258…; raw summary925/925 PASS, exit0, wall44.636s, sampled peak
process-tree RSS923983872B, sampled owned TMP0B, raw94B complete. Guard,
stop reason, overrun отсутствуют;592 completed samples. RSS/TMP sampling
не доказывает hard quota или transient peak и не покрывает произвольные
пути вне owned TMPDIR. Receipt SHA256:
`74b4c571b0d6cd9239cb26e3ef697a4b0a2d3164dd51f6db120d115dd250e49c`;
original raw log SHA256:
`918120b36de0f3b3d502d4cc5ba0e21497ba4b8f6ae32e6a895f48de6a35d7aa`.

Это наблюдение targeted manufactured functional run. Source CPP перед
execution был PASS, но independent **final CP-CODE и docs CPP pending**.
Healthy report обязан иметь только `{r_sum}` failure, N1 —
`{r_D,r_D_candidate,r_U,r_sum}`, N2 — `{r_D,r_D_candidate,r_sum}`.
Все эти expected report failures являются passing assertions fixture;
scientific acceptance, SCBA physical run, sufficient discretization,
experiment, S1 и fullCI/native acceptance из925 assertions не следуют.
GREEN-1 остаётся46pass/0fail/1error до coupled audits. Current counter2/2,
автоматического третьего D round или runtime rerun нет.

### Narrow D1 → D2 math mapping (frozen accepted artifact)

Следующий текст взят побайтно из approved mathdiff91959685… после его
верхнего заголовка. Его pre-execution формулировки — историческая стадия
design; текущие source/runtime statuses приведены выше.

Base D1 SHA256 `96c2c50a104b984c72603fa7066986d9bf6c02a8e59ba1ce0df9f06e321da162`.
Причина: GREEN-1 source guard exact wZ1 false, ни одного coupled audit.
D2 не меняет public builder arguments, source operators или S1; предполагается
только независимый D2 CP-DESIGN. Ни docs, ни test этой запиской не изменены.

| Привязка | D1 unit-spacing formula | D2 actual-spacing formula |
| --- | --- | --- |
| Expected spacing | h=1 без exact unit-chain proof | `h=RN(RN(RN(30*RN(1/1e9))/RN(10*RN(1/1e9)))/3)`; source prediction nextfloat1 |
| Fatal grid guard | `wZ==ones3` | exact `wZ==fill(h_ref,3)` plus exact canonical ell/L0/x chain; no isapprox fallback |
| x,wZ | `[.5,1.5,2.5]`,ones3 | `h*[.5,1.5,2.5]`,fill(h,3) |
| Input/grid length | Lp3 and grid3 treated identical | input Lp64 and real3h separate; no exact sum-period equality; r_n denominator uses input Lp64 |
| Discrete basis Φ | φ | φ unchanged (global phase only) |
| Continuous basis χ | φ | `φ/sqrt(h)` |
| Floating χ normalization | unit-spacing identity | source RN(sqrt(h64)) may1; actual Φ/χ and h-weighted normalization have explicit roundoff budget |
| Centre/Z,spread |1.5,1/√2 | `1.5h`, `h/√2` |
| FD hopping/H0 | t0,(2−√2)t0 | `tZ=t0/h²`, `(2−√2)tZ` |
| Embedding τ,c | −t0/4,t0²/8 | `−tZ/4`, `tZ²/8` |
| k kinetic | t0k² | t0k² **unchanged**, never tZk² |
| Donor density | D/3 | D/(3h) |
| Electron volume density | Dφ² | Dφ²/h |
| Sheet/charge | D | D unchanged; spatial integral uses h |
| Poisson ideal L,b | M,λD d/12 | M/h²,λD d/(12h) |
| Poisson floating L | M exact assumed | qL_ref M, `qL_ref=RN(1/RN(h_ref²))`, exact source-derived coefficient |
| U*,u* | λD d/36,λD/72 | λD h d/36,λD h/72 |
| Actual-coefficient Up reference | la D d/36 | la/(3qL_ref) P n_float |
| ζ actual reference | la(D−3nd)/(4√3) | la(sum(n_float)−3nd)/(4qL_ref√3) |
| ζ solve enclosure | short scalar rounding budget | h/q-corrected sameγ32 regression plus independent Q4 residual forward bound≤that budget, inverse2-norm4/3; no loosening |
| Neutrality reference | abs(D−3nd)/(3nd) | abs(sum(n_float)−3nd)/(3nd); source h-dot roundoff retained |
| Forward coefficient | 1/3 | 1/(3qL_ref) |
| Density SI conversion | nbar/L0³ | same units, but nbar=Dφ²/h; canonical L064 with conversion rounding |
| Roots z | E−H0−u−t0k² | E−H0_h−u_h−t0k² |
| Raw/normalized charge measure | gs/(2π)wEwk | **unchanged**, no h inserted |
| rK/A/Σ normalizers | full-array Frobenius+1e-14 | **unchanged** |
| rD normalizer | block `|Dg−1|/(|D||g|+1)` | **unchanged** |
| rP normalizer | norm(b)+la1e-14 | **unchanged** with h-correct b |
| Wrong-U rU | λD/18 | λD h/18 |
| Separation x | λD√8/(72t0) | λD√8 h³/(72t0) |
| rsum | max(abs(A/2π−1),1) | **unchanged** functional formula, h-correct A |
| Expected failed sets | H{rsum};N1{D,Dc,U,sum};N2{D,Dc,sum} | **unchanged** exact keys in full D2 |

`d=[−1,2,−1]`, `M=[2 −1 −1;−1 2 −1;−1 −1 2]`,
`P=I−ones3x3/3`, D1 constants/energy/radial grid remain fixed.

Future I delta is confined to: static h-dependent expected fields;
`independent_literals`; tZ vs t0 uses in `literal_state`/`hamiltonian_budget`;
actual-stencil/rounded-n Poisson enclosures; h-aware density/basis checks;
the two exact grid guards and donor/density expectations. Histories keep D
as sheet charge but carry h-correct density/U and same-field J. No new
production/formula/tolerance switch, optional N3, retry or scope extension.

γ16 complex products/γ32 six-entry sums/γ64 scalar Dyson remain unchanged;
the new h factors are exact input-parameter arithmetic evaluated outward,
with explicit short normalization/stencil formation bounds. Do not widen
these budgets to hide a wrong h dependency. D1 original artifact remains
failed evidence; D2 is a new candidate, never a rewrite of GREEN-1.

### Inherited exact API source bindings и ABI

Следующие public seams/source locators перенесены из accepted D1 без
изменений: D2 меняет h-dependent literals, а не подписи и record types.

* `src/numerics/reference/final_audit.jl:38`:
  `fresh_fixed_point_metrics(problem::NEGFProblem,scba::SCBAResult,
  U::AbstractVector{<:Real},candidate::FixedPointAuditCandidate)`.
  Публично экспортирован `fresh_fixed_point_metrics`; вводить новый API не нужно.
* `src/numerics/reference/validation.jl:516`:
  `validate_solution(solution::NEGFSolution;candidate=nothing,audit_failure=nothing)`.
* Там же `:751`: `QCLNEGF._validate_solution_metrics(solution::NEGFSolution,
  metrics;accept_approximate_scba=false)`. Это existing package-visible
  underscore seam, не новый stable public API.
* Там же `:33/:47/:60/:73`: определения rD/rA/rK/rΣ; `:636` fresh Poisson,
  `:688` limits. `r_D_candidate` имеет limit `t.r_K`, не `t.r_D`.
* `src/numerics/reference/basis.jl:14–84`: checked `:none` basis;
  `:288–318`: `H0+Φ†diag(U)Φ+Cbar*k²*Mminus1`.
* `src/numerics/reference/periodicity.jl:109–163`: literal target map
  `Σplus=Tplus*(Wplus G)*Tplus†`, minus аналогично; это **legacy bulk
  feedback**, не cavity embedding.
* `src/numerics/reference/poisson.jl:11–76`: cyclic stencil и bordered solve.
* `src/domain/types.jl:684/:777/:978/:1022/:1048/:1084`: complete records.

Pure record field inventory из frozen I2 source assertions (без численных
правил и новых defaults). Complete NEGFProblem fields:
`physical,numerical,scattering,scales,grids,profiles,basis,kernels,W₊ᴱᵖ,W₋ᴱᵖ,W₊ᴸᴼ,W₋ᴸᴼ,energy_shift_discretization,models`.

| Type | Frozen field names |
| --- | --- |
| `GreenState` | `Gᴿ,Gˡ,Gᵍ,A,condition_number,dyson_scale` |
| `SelfEnergyFamily` | `Σᴿ,Σˡ,Σᵍ` |
| `FixedPointAuditCandidate` | `scattering,embedding,plus,minus,roundoff` |
| `SCBAResult` | `green,scattering,embedding,embedding_plus,embedding_minus,history,converged,status,quality,restart_contract,mixer_state` |
| `NEGFSolution` | `problem,options,Uᴴ,n,scba,outer_history,observables,report,converged,status` |
| `OuterIteration` | `μ,r_P,r_U,r_n,r_neutral,r_J,r_Jchange,r_population,ζ,J,populations,density` |
| `SCBAIteration` | `ν,r_D,r_A,r_K,r_Σ,r_λ,λ,r_PSD,r_caus,r_roundoff,r_Jchange,r_population,J,raw_charge,target_charge,normalized_charge,lambda_change,witness,physical_markers` |

## O08-04 — nonzero public current и coherence

**API/units.** `energy_resolved_current(problem::NEGFProblem,
family::SelfEnergyFamily,green::GreenState)`;
`boundary_flux(problem,family,green)`;
`bond_current(problem::NEGFProblem,Gl::AbstractArray{<:Number,4})` —
`src/physics/observables.jl:224,240,270`. Green/self-energy arrays scaled,
shape `(NE,Nk,Nb,Nb)`. Returned jE: A m⁻² eV⁻¹; flux: m⁻² s⁻¹;
internal bond currents: A m⁻², Nz−1 entries, seam отсутствует.
Outgoing engineering electron-flow current `J=e Φ`, conventional charge
current имеет противоположный знак.

**Complete public carrier.** Checked helper выше с `Nz3,Nb2,:none`,
`NE3,Nk2,P1,Nφ4,Nqz3`; actual `ε=[−1,0,1]`, `wE=[1/2,1,1/2]`,
`κ=[0,1]`, `wk=[1,3]/(16π)`, gs2. Ранее предложенные NE2/Nk1 и
свободные dyadic weights заменены: public constructor их не допускает.
Трёхузельная uniform Dirichlet цепочка даёт independent basis

```text
Hcell=t*[2 -1 0;-1 2 -1;0 -1 2],
φ1=(1/2,1/√2,1/2), E1=(2−√2)t,
φ2=(1/√2,0,−1/√2), E2=2t, Φ0=[φ1 φ2], χ0=Φ0.
```

Не фиксировать tied-endpoint phase LAPACK для φ2. После independent
assertion `abs(dot(φa,Φactual[:,a]))≈1`, определить
`d[a]=dot(φa,Φactual[:,a])/abs(dot(...))`, P=Diagonal(d), проверить
`Φactual≈Φ0P`, `χactual≈Φ0P`. **Все** literal Green и family matrices
перевести как `Xactual=P'XliteralP`. Actual Φ используется лишь для
column-phase transport; expected current не вычисляется из actual Φ/Hgrid.

Полный remaining basis inventory (с phase transform P при необходимости):

| Field | Independent value / объяснение |
| --- | --- |
| localization | `:none` |
| centres,spreads | `[3/2,3/2]`, `[1/√2,1]` |
| overlap eigenvalues,κ | ones2,1 within arithmetic budget |
| H0 | diag(E1,E2) |
| Z | `P'*[3/2 -1/√2;-1/√2 3/2]*P` |
| M−1 | I2/mr |
| T+,T− | `−t*P'*[1/4 1/(2√2);−1/(2√2) −1/2]*P`, adjoint |
| window_eigenvalues,Egrid0 | `[E1,E2]`,2t; :none stores central-cell energies |
| r_orth,r_eigen | builder-computed, finite and existing static gates; never overwrite0 |
| three r_translation fields | builder0 for :none, outside PzP translation gates |

**Nonzero literal fixture.** Все arrays shape `(3,2,2,2)`. В analytic basis
в каждом (e,m):

```text
Q=[2 i;−i 1], B=[1 i;−i 2], F=I2/2, D=I2,
a=[1 2;3 4;5 6] indexed(e,m),
G<=iaQ, G>=−iI2, GR=−i(aQ+I2)/2, A=aQ+I2,
Σ<=iI2/2, Σ>=−iB, ΣR=−i(B+I2/2)/2.
```

Все correlation matrices PSD, обе Keldysh identities заданы явно. Это
current functional input, не Dyson/SCBA solution. Independent scalar sums:

```text
tr(B*aQ−F*I)=a*(2+1+1+2)−(1/2+1/2)=6a−1,
S=[5 11;17 23;29 35],
sum_m wk*S=[38,86,134]/(16π)=[19,43,67]/(8π),
jE=[19,43,67]/(8π²) * J0ref/E0eV,
sum_e (E0eV*wE[e])*jE[e]=43/(4π²) * J0ref,
eSI*boundary_flux=43/(4π²) * J0ref.
```

J0ref/t используют independent SI literals из carrier section, не actual
scale readback или CODATA. Сначала отдельно проверить scales against SI,
затем оба public API против retained expected values и Unitful units.
Physical dE=`[.05,.1,.05] eV`, не dimensionless wE. Общий ошибочный
multiplier/sign/2π двух APIs не сокращается в одном agreement.

**Complete Green diagnostics.** GR здесь намеренно не решение заданного
Dyson D; family — externally supplied manufactured current input, не actual
scattering/embedding map. `condition_number`/`dyson_scale` всё же относятся
к declared actual D, без unrelated ones/zeros. В analytic basis
`hjm=Ej+t*κm²`, `ξj=εe−hjm`:

```text
D_em=[ξ1+3i/4 −1/2;+1/2 ξ2+5i/4],
τ=ξ1²+ξ2²+21/8,
d2=(ξ1ξ2−11/16)²+((5ξ1+3ξ2)/4)²,
q=sqrt(max(0,τ²−4d2)),
condition_number[e,m]=sqrt((τ+q)/(τ−q)),
dyson_scale[e,m]=max(sqrt(ξ1²+9/16)+1/2,sqrt(ξ2²+25/16)+1/2).
```

Imaginary Hermitian part D=(B+F)/2 positive definite, D nonsingular.
Discriminant below0 beyond inherited arithmetic budget → stop, не PSD
repair; max0 только для roundoff отрицательного discriminant.
До required plain `GreenState`/`SelfEnergyFamily` construction test-local
checking assembly проверяет shape/ComplexF64/finite всех четырёх Green и
трёх family arrays, metadata `(3,2)` finite/positive по формулам, no-alias,
phase mapping, A Hermiticity, обе Keldysh identities, PSD signs.
Это all-field functional check, не сертификат Dyson/SCBA solution.
Не получать expected arrays через Green/self-energy helpers.

**Coordinate reconstruction — те же actual axes.** Для G<=iaQ и Φ0:

```text
Re[(Φ0*(iQ)*Φ0')[2,1]]=Re[(Φ0*(iQ)*Φ0')[3,2]]=−1/2,
sum_em wE[e]wk[m]a[e,m]
 =((1+6)/2+(3+12)+(5+18)/2)/(16π)=15/(8π),
H12=H23=−t,
expected bonds=[1,1] * 15t/(8π²) * J0ref.
```

Равенство этих двух bonds выведено только для данной sine basis/Gl;
universal uniform local J не требуется. Periodic seam отсутствует.
Diagonal-only Gl сохраняет G11=2ia,G22=ia и даёт zero bonds.
Для truncation отдельный checked builder `(3,1,:none)` с теми же physical
axes; проверить sole column φ1 up to phase, подать `(3,2,1,1)` Gl=2ia.
Его reconstructed entries purely imaginary, bonds0. Не оставлять Nb2
metadata, не мутировать basis, не перенормировать density и не называть
этот truncated carrier той же физической моделью.

**Норма/допуск/controls.** Existing finite-chain current regression budget
`atol=2e-14,rtol=2e-13` применяется **после деления actual/expected current
на independent J0ref**, resolved — J0ref/E0eV. Такой же small arithmetic
budget применяется к dimensionless analytic basis entries/scale ratios;
existing static gates сохраняются. Это regression policy, не2e-14 A/m²
и не LAPACK-specific forward-error theorem. Near zero absolute branch.
Expected scalar sums dyadic exact до π/scale, radicals/SI conversion
проверяются тем же budget. Sign/omitted2π/dropwk/dropwE/coherence/truncation
mutations обязаны отличаться от retained expected values:
dropping wk1 loses `[5,17,29]/(16π²)` from scaled jE; dropping wk2 loses
`3[11,23,35]/(16π²)`. Omitting wE gives integral129/(8π²), versus86/(8π²),
ratio3/2; dimensionless weights as physical dE give ratio10. Dropping family
coherences changes trace6a−1 to4a−1. RHS не менять вместе с mutation.

**Primary/reuse.** Wacker `0706.4202v1`, Eqs. (3),(5),(6), Fig.1 объясняют
coherence/finite-basis caveat; независимые numerical RHS выше выведены
непосредственно из stated finite matrices. Existing cavity/Landauer test
проверяет finite_chain_green, не public bond current; повторять его не нужно.

## O08-05 — absolute scattering scale: design_blocked

**API/shape.** `lo_phonon_kernel(p::PhysicalParameters,s::ScaleSystem,
grids::ModelGrids,basis::BasisData) -> (K,FLO)`;
`interface_roughness_kernel(p,s,grids,basis) -> K`;
`acoustic_kernel(p,s,grids,basis) -> K`;
`alloy_kernel(p,s,grids,profiles::MaterialProfiles,basis) -> K`.
`src/numerics/reference/kernels.jl:201,240,280,324`.
`K` shape `(Nk,Nk,Nb,Nb,Nb,Nb)` with indices `(m,mprime,a,c,d,b)`;
raw scaled K, не normalized `KernelSet.K`; FLO shape `(Nq,Nb,Nb)`.
`K_SI` units J²m², divide by independently computed `(E0_J*L0_m)²`.
χ dimensionless, ψ=χ/sqrt(L0); wˣ dimensionless, qbar=q*L0;
angular kernel average has 1/Nφ, radial wk лишь в contraction.

| Channel | Owning RHS и literal geometry candidate | Primary gap / status |
| --- | --- | --- |
| LO | `e²*hbaromega/(2ε0)*(1/epsinf-1/epss)` [J²m] times `∫dqz/(2π) FF*/(q²+qz²+qs²)` [m]; one-state x=0,wˣ=1,χ=1; qbar=2,qs=0; qzbar=(-1,0,1),wqbar=(1/2,1,1/2), F=1; integral `9L0/(40π)` | Lee–Wacker v1 B20/B23 establish unscreened prefactor with **absolute** ε; relative ε in Core introduces ε0 explicitly. Arithmetic candidate derived; one-point envelope is a discrete coefficient fixture, not continuum geometry validation. Complete parameter carrier and roundoff budget still pending |
| IFR | `πΔ²Λ² exp(-q²Λ²/4) * (ΔEc)² * |ψ(zj)|⁴` for one independent interface; Δ/Λ in m, ΔEc in J, ψ² in m⁻¹ | Primary Lee–Wacker B8 uses **exponential** correlation; it cannot certify Core Gaussian spectrum by substitution. Need explicit Gaussian covariance Fourier derivation and real periodic interface list/seam once/χ interpolation fixture; no default band jump guesses |
| acoustic | `Xi²*kB*T/(rho*vs²) * ∫|ψ|⁴ dz`, normalized constant envelope length ℓ gives integral1/ℓ | Author review 1412.3563 Eqs166–171 gives elastic/equipartition derivation, but primary-source convention/normalization chain not completed here. Lee–Wacker B26 uses fitted strength and is not absolute oracle |
| alloy | `Omega0*(DeltaV)² * ∫x(1-x)|ψ|⁴ dz`; constant x and length ℓ gives x(1-x)/ℓ | Need primary Bernoulli-site covariance/volume convention and mapping of DeltaV to stated material contrast; review204–205 alone not accepted primary proof here |

Owning equations `docs/src/theory/08_kernels.md` EQ-KLO/KIFR/KAC/KALLOY
дают точную модель для сравнения, но копия target formula не независимое
подтверждение её prefactor. Знать SI dimensions недостаточно для numerical
factor 2 или 2π. LO local algebra ещё не разрешает DEV07 всех каналов.

**Нормы/near-zero/tolerance gap.** Healthy absolute magnitude comparison
должна быть componentwise K_SI, потом residual с independent nonzero Kref;
при exact zero use absolute error, не делить на Kref. Требуется explicit
roundoff/quadrature bound для принятой geometry. Нынешние qK*K roundtrip
(`test/physics/exact_kernel_rescaling.jl`) и PSD tolerances не являются
absolute physical tolerance. Такой derived budget для всех четырёх сейчас
не получен; нельзя произвольно выбрать rtol для нового physical oracle.

**Positive/negative frozen intent.** Healthy positive prefactors при
positive measures; global K→2K обязан нарушать retained absolute RHS,
даже если normalized LO KMS/balance остаётся. IFR Δ=0 или ΔEc=0 →exact0;
acoustic T=0 →algebraic0 (не проверка применимости equipartition при0K);
alloy x=0/1 или missing required params →exact0. LO missing 1/(2π) и
changed ε convention — обязательные negatives. Unknown primary inputs не
подставлять нулём. Existing impurity midpoint fixture сохраняется отдельным
каналом и не закрывает четыре gaps.

**Next design.** Разделить четыре channel passports на finite child
packets с primary formula+version+hash, permitted geometry и arithmetic
bound. До их CP-DESIGN DEV07 не READY; model/threshold change запрещены.

## O08-06 — PV: discrete/direct/FFT/continuum раздельно

**Exact APIs.** `direct_hilbert_transform(Γ::AbstractArray{<:Number,4},
ε::AbstractVector{<:Real},wE::AbstractVector{<:Real};return_roundoff=false)`
— `reference/kernels.jl:500`; `fft_hilbert_transform(Γ,ε,wE;
column_chunk::Integer=32,return_roundoff=false)` с теми же type annotations
— `optimized/production_operators.jl:723`. Shape Γ `(NE,Nk,Nb,Nb)`, Λ
того же shape/units; ε,wE в E0 и мера/denominator сокращают E0.
Никаких k/z weights. FFT требует uniform increasing ε, positive chunks,
linear zero padding L=nextpow(2,2NE-1), extraction noncyclic. Direct map
`Λ_i=Σ(j≠i) w_j Γ_j/[2π(ε_i-ε_j)]`, диагональная cell опущена.

**Independent discrete RHS.** NE=3,Nk=Nb=1,

```text
ε=[-1,0,1], w=[1/4,1/2,3/4], Γ=[1,2,4]
expected Λ=[-5/(4π), -11/(8π), 9/(16π)].
```

Каждую сумму раскрыть ровно в две scalar terms; expected literal array не
вызывает direct/FFT/kernel builder. Оба APIs отдельно сравнить с этим RHS,
не только друг с другом. Shapes/weights не physical quadrature campaign,
а discrete-map sentinel. Wrong denominator sign, finite nonzero diagonal
addition или cyclic-wrap term меняют retained nonzero endpoint values.
`Γ=0` supplemental control даёт Λ=0; diagonal contribution mutation с
division by0 обнаруживается finiteness + all-entry comparisons.

**Отдельный finite-support continuum oracle.** Manufactured scalar
`Γ(x)=4+x` на `[-3,3]`, ноль снаружи. В target E=0 positive Γ на всём
носителе; singular constant term cancels exactly in PV, linear term gives
`PV∫[-W,W] Γ(x)/(2π(0-x)) dx = -min(W,3)/π` для W=1,3,5.
Полный infinite-domain value с declared compact support равен `-3/π`.
Tail omitted for W=1 is `-2/π`, for W≥3 exact0; это точный tail для данного
manufactured function, не bound физического SCBA spectrum.

For odd N=17,33,65 use `h=2/N`, window W=1,3,5 with respectively
N,3N,5N cell-centered nodes, all weights=h. Target 0 — middle node;
support boundaries ±3 — cell edges. На target каждая symmetric constant
pair cancels; каждое noncentral supported cell даёт ровно `-h/(2π)`.
Отсюда независимо от target helper

```text
Λ_h(0) = -min(W,3)/π + h/(2π)
continuum quadrature error = h/(2π), positive, ratio h1/h2
full-support error = Λ_h(0) - (-3/π)
                  = (3-min(W,3))/π + h/(2π).
```

Это exact center-only discretization error/bound при stated quadrature;
ошибку окна и опущенной cell можно asserted отдельно. При fixed h
W1→W3 измеряет known tail; W3→W5 не меняет support и ловит wrap/
padding contamination. При fixed W refinement N даёт exactly derived
trend; не требовать second order от omitted-cell operator. После вычитания
known analytic h/(2π) остаток ограничивается arithmetic regression budget.
General all-energy continuum PV bound **не выведен** и не заявляется.

**Норма/tolerance/near-zero.** Componentwise scaled Λ comparison с existing
`test/numerics/fft_linear_hilbert_map_equals_direct_principal_value.jl`
`rtol=atol=2e-13` для discrete literals и center continuum formula;
absolute branch нужна zero control. Это inherited small-array FFT/direct
regression policy, не continuum error threshold. По малым discrete sums
можно дополнительно дать γ_n absolute bound с explicit operation count;
он не переносится автоматически на FFTW. Continuum error по exact formula
не «проходит» из r_roundoff. Returned r_roundoff direct reverse-order и
FFT doubled-padding имеет норму `||ΔΛ||F/(||Λ||F+1e-14)` без quadrature;
это sensitivity diagnostic, не требуемая continuum bound.

**Accepted clarification — propagated ratio budget, CP-DESIGN/CP-CODE PASS.**
Targeted DEV08 GREEN-1 (fixture SHA256
`bae8a57957cd9df4da0d7369521a7d7ffa776b4c095859e0d1a9f86a9607eb6e`)
дал13/13 discrete и107 pass/1 fail continuum: ratio65/33 measured
`1.969696969698286` против `1.9696969696969697` при2e-13. Original failed
log SHA256 `495227f26beba7ca0b274de6e8f39bf6ada9e1015afbe21e41a7f000fb250fa7`,
receipt SHA256 `e0f7a3b5a8b2afb987232333d04cb3ab6a2600de25dc871f388af5836b8ef775`.
Сохранить original failed evidence в root execution/evidence, не заменять
успешной повторной попыткой. Этот факт не доказывает defect PV operator.
Current fixture4533b105… после accepted correction имеет GREEN2
211/211 assertions (13 literal,198 continuum/scalar), matching raw
receipt/log identities перечислены выше. Historical GREEN1 остаётся failure.

Принятая correction меняет только propagation для quotient, **все
underlying functional comparisons и ATOL=RTOL=2e-13 сохраняются**.
Julia1.13 scalar isapprox (`base/floatfuncs.jl:222–228`) использует
`abs(x−y)≤max(ATOL,RTOL*max(abs(x),abs(y)))`, не сумму budgets.
Для `0<x,y<1`, ATOL=RTOL, это абсолютная ветвь A=2e-13.

Пусть independently `a=eA=(2/NA)/(2π)`, `b=eB=(2/NB)/(2π)`,
measured x,y проверены отдельно этим же scalar gate. Их отношение
имеет conditioning

```text
| (a+δa)/(b+δb) − a/b | ≤ (B+(a/b)B)/(b−B),
|δa|,|δb|≤B, b>B.
```

`B=nextfloat(A)` включает rounding floating subtraction в isapprox,
не является новым PV tolerance. При e65≈.0049 деление усиливает absolute
budget до порядка1e-10; прямой ratio rtol2e-13 не следует из functional
policy. Для actual predicate использовать outward bounds, без empirically
chosen tolerance:

```text
a_lo=prevfloat(a−B), a_hi=nextfloat(a+B),
b_lo=prevfloat(b−B), b_hi=nextfloat(b+B),
r_lo=prevfloat(a_lo/b_hi), r_hi=nextfloat(a_hi/b_lo).
```

Каждый endpoint/result finite; `0<a_lo≤a_hi<1`, `0<b_lo≤b_hi<1`.
Если denominator lower≤0, stop unresolved ratio; floor не добавлять.
Outward addition/subtraction/division и monotonic positive Float64 division
enclose computed x/y. Independent mathematical ratio ρ=NB/NA также должен
лежать в [r_lo,r_hi]; π/h literals rounding не скрываются в expected actual.

Accepted frozen predicates ниже — псевдокод contract; reviewed implementation
реализует также требуемые explicit fatal guards перед dependent division:

```julia
nodes = (17,33,65)
expected_errors = [(2/N)/(2π) for N in nodes] # independent, not helper output
for i in eachindex(errors)
    @test isfinite(errors[i]) && 0 < errors[i] < 1
    @test 0 < expected_errors[i] < 1
    @test isapprox(errors[i],expected_errors[i];atol=ATOL,rtol=RTOL)
end
@test errors[1] > errors[2] > errors[3] > 0 # retained
for (i,j) in ((1,2),(2,3))
    a,b = expected_errors[i],expected_errors[j]
    B = nextfloat(ATOL)
    alo,ahi = prevfloat(a-B),nextfloat(a+B)
    blo,bhi = prevfloat(b-B),nextfloat(b+B)
    # failed finite/positive precondition must stop before division
    @test all(isfinite,(alo,ahi,blo,bhi))
    @test 0 < alo <= ahi < 1 && 0 < blo <= bhi < 1
    @test bhi < alo
    lo,hi = prevfloat(alo/bhi),nextfloat(ahi/blo)
    @test isfinite(lo) && isfinite(hi) && 0 < lo <= hi
    @test lo <= nodes[j]/nodes[i] <= hi
    ratio = errors[i]/errors[j]
    @test isfinite(ratio) && lo <= ratio <= hi
end
```

`@test` не short-circuit: reviewed implementation останавливает dependent
division при failed precondition, не продолжает с invalid endpoints. Три per-error
assertions относятся к отдельно пересчитанным error-loop values, а не
предполагают equality с более ранним W3 operator invocation.

Cheap scalar controls без новых operator evaluations: a,b и допустимые
perturbations `a+A/2,b−A/2` проходят strict functional gates и interval;
constant ratio1 и second-order ratioρ² лежат вне first-order interval.
Common factor2 или общий sign reversal сохраняют ratioρ: их отвергают
отдельные amplitude/positivity gates, а ratio сам их не диагностирует.
Явно проверить эти accept/reject утверждения при retained RHS. Ratio —
consistency evidence known omitted-cell refinement, не независимый
сертификат physical discretization или substitute signed Λ oracle.

После independent high CP-DESIGN medium заменил только две ratio assertions
и добавил эти same-budget guards/controls; CP-CODE PASS подтверждён receipt
2cc83ef3…. Captured GREEN-1 bytes сохраняются; GREEN2 потребил оставшуюся
invocation existing DEV08 L-JL, не новый счётчик. Underlying
Λ/physical tolerances, operators, windows, meshes и FFT policy не менять.

**Reuse/граница.** Existing direct Gaussian odd-symmetry и FFT/direct
random-array anchors остаются. `analytic_product_integration.jl` и frozen
triangle покрывают piecewise-linear integration с exterior taper, иной
оператор и boundary model. Здесь им не приписывается direct/FFT pass.

## Primary-source receipt и стоимость

Три PDFs прочитаны через web (arXiv versions видны в тексте), скачаны
в /tmp только для hash/provenance; aggregate 3,682,827 bytes. Каждый curl
имел max-filesize32 MiB, connect timeout10s, total timeout30s; общий paper
ceiling96 MiB≤128 MiB. Первая sandbox DNS попытка завершилась exit6;
bounded read-only network escalation с тем же первичным URL дала exit0.
Никакой network installation. Сохранить URLs/hashes в owning document,
не коммитить PDFs. Source formulas перепроверить по original page если
web math extraction неоднозначна.

| Source URL/version | SHA256 downloaded PDF | Formula scope |
| --- | --- | --- |
| https://arxiv.org/pdf/cond-mat/0212059v1 | dfabd791bfb2aa4603233e9b545879e07a3bfc63e8906cbc810760726358ea0c | Keldysh relation; IFR vertex B10/B11 with exponential covariance B8; bulk unscreened LO B20/B23 |
| https://arxiv.org/pdf/0706.4202v1 | 091a39824c3fb4b5f1f3b0f0659d4a99663b48b6949806aa8586710b6d1d937b | local current Eqs3,5,6 and Fig1 finite-basis caveat |
| https://arxiv.org/pdf/1412.3563v1 | e2e43e27c96994bfa5aa8750044bef33afe92de1360c0c0c5861bf235d8ff7c8 | secondary author review; source tracing only, acoustic166–171, IFR197–199, alloy203–205; not promoted to primary pass |

An unrelated arXiv URL `1302.3725` was accidentally opened during browsing,
identified as spin-ice paper and excluded immediately; no equation/evidence
from it is used. `cond-mat/0109163v1` was read as historical context only;
its phonon text extraction is not used as absolute kernel oracle.

## Завершение source/design checkpoint

| Утверждение | Свидетельство | Ограничение | Решение |
| --- | --- | --- | --- |
| Keldysh/current/PV получают различающие RHS | Literal matrix/component sums, analytic PV equations, accepted design/code и658 targeted assertions | Bounded manufactured functional scope, не coupled science | сохранить |
| Poisson incompatibility видна через ζ/rP/rneutral | Constant RHS bordered equations, exact-ND precondition, existing sine derivative anchor | Assessment carrier synthetic, full coupled state не доказан | сохранить |
| Complete functional final U/G D2 принят | D2eb9e1a8c…/review53c3e9b4…; I2 81fbb258… source CPPbb5b5260… PASS; executor GREEN-2 925/925 targeted assertions | Independent final CP-CODE/docs CPP pending; healthy global report FAIL r_sum; S1/fullCI/native не приняты; GREEN-1 сохранён | дополнительно измерить |
| Four absolute prefactors сертифицированы | Owning units/formulas и partial primary LO chain | Geometry/tolerance/primary gaps, DEV07 blocked | данных недостаточно |
| Current PV targeted runtime прошёл | GREEN2 fixture4533b105…,13+198 assertions; historical GREEN1 failure сохранён | FullCI/SCBA/science0 | сохранить |
| SCBA/science выполнены | Source/design, paper hashes, четыре bounded functional runs | Convergence/discretization/experiment не установлены | данных недостаточно |

Независимый CP-DESIGN принял пять bounded functional scopes;
O08-05 остаётся design_blocked. Reviewer проверил literal arithmetic,
shape/units, applicability, fixed RHS mutations и статусы. Автор документа
не закрывает собственный CP-CODE. Отдельный I medium выполняет только
accepted scope с envelope соответствующего DEV03/04/05/06/08; PV ratio
predicate получил CP-DESIGN/CP-CODE и bounded GREEN2 evidence. DEV05 I2 source PASS;
GREEN-2 925/925 functional assertions PASS, independent final CP-CODE/docs
CPP pending. DEV07 не READY.
Полная six-oracle readiness этим документом не устанавливается.

### DEV07 — accepted D1 absolute covariance fixture (2026-10-19)

Этот адресный subsection заменяет `design_blocked` **только для описанной
ниже новой fixture**; исторические candidates выше не являются её inputs.
`test/physics/absolute_scattering_kernel_oracles.jl` обнаруживается existing
physics group автоматически. Source authoring не является runtime receipt.
Production kernels, S1 model и scientific acceptance thresholds не меняются.

Принятый D1: SHA256
`1d74158042e305c37875ceef25589b08f6ca9256fdd7f65993ee6825c165a07e`;
независимый CP-DESIGN: SHA256
`caea42e27026e6bee60d1f59988a4f1a927820c503f50b2c07f8a0803ba52687`.
Owning source — Core `7324d98c0fccd61b4a1a2756b5cabe2fd0e0e49c`.

Checked builders создают три physical cells длины `ell=1e-8 m`, массы `.067`,
band edges `[2,5/2,1/2]tau`, где `tau=hbar²/(2m0*.067*ell²)`, и Nb=1,
localization=:none. Literal central block в units tau равен
`[4 -1 0;-1 9/2 -1;0 -1 5/2]`, eigenvector `[1,2,4]/sqrt(21)`, eigenvalue2,
gap3/2. Независимый interval residual и norm certificate обязательны;
результат eigensolver не подставляется в expected amplitude.

`k*ell=[0,2]`, phi=[0,pi/2,pi,3pi/2], qz*ell=[-pi,0,pi].
Для любого f(v), v=q²ell², full radial2×2 RHS равен
`[f(0) f(4);f(4) (f(0)+2f(8)+f(16))/4]`.
Все четыре положительных raw entries имеют SI units J²m²; каждая scale
`.1eV/10nm` и `.2eV/20nm` отдельно делит literal RHS на `(E0_J*L0_m)²`.
Ни actual F, ни CODATA/scales readback, ни normalized K не строят expected.

| Канал | Независимый discrete SI f(v) | Primary normalization и граница |
| --- | --- | --- |
| LO | `e²*(.05e)/(32eps0) * ell/2 * [1/(v+1)+(169/441)/(v+1+pi²)]`; literal F=`[-13i/21,1,13i/21]` | [Lee–Wacker v1](https://arxiv.org/abs/cond-mat/0212059v1), B19–B23: bulk Fröhlich/absolute permittivities/dqz/(2pi). Positive screening1/ell — owning EQ-KLO assumption; finite trapezoid, не infinite integral |
| IFR | `pi*(1e-10)²*tau²*(2329/3136)*exp(-v/4)` | [Huang et al. v2](https://arxiv.org/abs/2205.00365v2), Eq8 covariance convention; direct 2D Gaussian transform. Lee–Wacker B10–11 band-step vertex. Independent interfaces/seam once; midpoint D*ell/tau=`[1/14,3/56,-6/7]`. Carrier screening не добавляется |
| acoustic | `(4e)²*kB*70/(4000*4000²)*13/(21ell)` | [Archetti–Neophytou v1](https://arxiv.org/abs/2012.12344v1), Eq20/Appendix: elastic equipartition sum of emission+absorption cancels vertex1/2. Planar projection derived independently of nanowire geometry |
| alloy | `3.125e-29*(.6e)²*13/(112ell)` | Huang Eq9: Bernoulli covariance `Omega*DeltaV²*x(1-x)`; Omega=a³/4 for a=.5nm, x=1/4. Synthetic contrast is independent of band jumps; bare disorder model |

Frozen primary PDF SHA256 respectively:
`dfabd791bfb2aa4603233e9b545879e07a3bfc63e8906cbc810760726358ea0c`,
`00297abb84fdc44fc865461b9ba102dfef7ff7cde2924b2999b0deca054f5f5c`,
`843184eec1632f8f46fde7e3a2492126763487dc99172f0950761239f9ae046e`.
Source authoring performs no new literature retrieval or device-model validation.

Reference evaluation uses 256-bit decimal-string intervals, with outward
rounding per primitive (including products, reciprocal, square, sqrt, exp).
`u=2^-53`, `g=gamma(32)`, `a=gamma(1024)` and the exact D1 §4 propagation
formulas define per-channel delta; no empirical epsilon is fitted.
Grid/profile bounds and independent residual/gap bound beta protect shared
carrier code. Scalar cos/sin/exp are certified against MPFR at each actual
represented Float64 argument, including the complex LO exponential outputs;
sqrt assumes correctly rounded IEEE normal arithmetic. Failed certificates
remain failures and may not be handled by widening delta.

Operation inventory counts real rounded operations (complex multiply6,
complex add2); libm calls have their own one-ulp certificates. Independent
CP-CODE supplies conservative ceilings per raw entry: LO744, IFR544,
acoustic268, alloy280, all <992<1024. The earlier I1 comment understated
the explicitly counted LO denominator operations; this E1 accounting
erratum changes no executable expression, gamma_n or tolerance.

| Channel | Independent real-op ceiling | Expansion |
| --- | --- | --- |
| LO | 744 | SI128 +9 envelope/phase terms×24 +4 angles×(q12 +3 qz blocks×24 +scale12) +accumulation16 |
| IFR | 544 | SI128 +3 interfaces×32 +4 angles×(q12 +spectrum16 +3 vertex contractions×12 +scale12) +accumulation16 |
| acoustic | 268 | SI128 +3 quartic terms×32 +prefactor16 +scale12 +extra accumulation16 |
| alloy | 280 | SI128 +3 variance/quartic terms×36 +prefactor16 +scale12 +extra accumulation16 |

SI128 includes canonical conversions/literal inputs and conditioned band
jumps (condition9). A LO qz denominator has three squares and two adds;
four factor operations and complex products/accumulation10 give19<24.
The accepted D1 universal ceiling992 and gamma1024 remain unchanged.
Grid/basis errors are separately certified. LAPACK is excluded from the
product count because its returned vector must pass the independent
residual/gap certificate, rather than inheriting an assumed LAPACK bound.

Each retained positive interval [Rlo,Rhi] checks real K in
`[Rlo*(1-delta),Rhi*(1+delta)]` and `abs(imag K)<=Rhi*delta`.
`delta<1/8` guarantees rejection of global2; normalize(2K)=normalize(K)
is only an algebra demonstrator. Per-channel prefactor/Gaussian/seam mutants,
positive screening2/ell and T35 sentinels, exact-zero IFR Delta0 and alloy
x0/x1/disabled cases are separate controls. No T0 or zero LO-screening
constructor is used. The impurity midpoint anchor and exact-kernel-rescaling
test remain unchanged; their narrower independence is not upgraded.

Runtime, SCBA/Poisson convergence, physical gates, continuum refinement and
experimental validation each require their own evidence. Nb=1 cannot test
general four-index orientation or adjoints; no such coverage is claimed.
