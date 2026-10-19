# DEV02 — паспорта независимых oracle O08-01…06

Статус: **CP-DESIGN PASS для bounded scopes O08-01/02/04/06**
(2026-10-19, независимый reviewer high). O08-03/O08-05 остаются
**design_blocked**. Полная готовность шести oracle и полное завершение
DEV02 acceptance **не установлены**.

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
| O08-03 | design_blocked | Не завершён complete manufactured final-state oracle с независимыми U/G/n/candidate/history |
| O08-04 | design_accepted | Actual NE3/Nk2 public carrier, signed SI current, literal coherence и genuine Nb1 truncation |
| O08-05 | design_blocked | Четыре absolute channel passports ещё не имеют законченной primary→geometry→SI→tolerance цепочки; LO derivation — локальный кандидат |
| O08-06 | functional/ratio design_accepted, CP-CODE PASS | Literal direct/FFT, center-only continuum/tail и accepted derived ratio interval; targeted GREEN2 passed |

`design_accepted` означает только прошедший независимый CP-DESIGN
bounded functional-regression scope. DEV03/04/06/08 имеют текущую
CP-CODE и targeted runtime evidence в своих envelopes. Runtime acceptance
не следует из design; ниже приведены exact fixture/raw receipt identities.
DEV03 ← O08-01; DEV04 ← O08-02; DEV05 ← O08-03
(blocked); DEV06 ← O08-04; DEV07 ← все четыре O08-05 (blocked, допустим
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

## O08-03 — complete final U/G: design_blocked

**Существующие seams.**
`fresh_fixed_point_metrics(problem::NEGFProblem,scba::SCBAResult,
U::AbstractVector{<:Real},candidate::FixedPointAuditCandidate)` —
`src/numerics/reference/final_audit.jl:34`;
`validate_solution(solution::NEGFSolution;candidate=nothing,audit_failure=nothing)`
— `validation.jl:516`, затем `_validate_solution_metrics`.
`project_hamiltonians(problem,U)` использует `Φ†diag(U)Φ`;
fresh rD/rD_candidate именно этот U. Final Poisson rU сравнивает freshly
solved Up с `solution.Uᴴ` (`validation.jl:636`).

**Неподменяемый контракт.** Требуется complete `NEGFProblem` с physical,
numerical, scattering, scales, grids, profiles, basis, kernels, four shift
maps, discretization, models. `GreenState` содержит GR/G</G>/A/condition/
Dyson scale. `SCBAResult` содержит все enabled scattering families,
embedding/+/- и SCBA history/classification; `FixedPointAuditCandidate` —
fresh same-state enabled families, embedding/+/- и measured roundoff. `NEGFSolution`
добавляет published U,n, outer history, options и report.
Все densities/charges/gates должны иметь объяснение из одного состояния;
unrelated zeros или arbitrary partial state не complete assembly oracle.

**Independent RHS/gap.** Нужны одновременно: literal Dyson/Keldysh family,
independent Poisson n(U), target raw charge, explicitly defined raw-map
candidate и histories с prescribed applicable gates. Scope manufactured
limit и его отличие от production baseline надо определить отдельно;
само требование alignment не доказывает, что все пять каналов обязательно
должны иметь ненулевой amplitude в test fixture. В данном source/design
пакете такая complete manufactured coupled state не выведена. Existing
`fresh_candidate_residual.jl` вызывает только `_dyson_residual`, а causal
quadratic root/maxiter harness — primitive SCBA components; ни один не
вызывает complete final assembly. `_reference_audit_candidate` использует
общие physical kernels и embedding, поэтому сам не independent RHS.

**Заданный будущий negative.** Для derived healthy complete state сохранить
G и изменить U на δU с `Φ†diag(δU)Φ≠0`, либо при fixed U оставить stale G;
retained independent Dyson denominator/Poisson candidate обязаны нарушить
fresh designated gate даже при stored small Σ residual. Нормы:
max-block `||DG-I||F/(||D||F||G||F+sqrt(Nb))`, без weights;
rK global Frobenius/(||G<||F+1e-14); rU absolute ∞ norm в E0.
Допуски — действующие rD/rK/rU, но quantitative positive/negative margin
ещё **не получен**, потому DEV05 не READY.

**Самый дешёвый next design.** Отдельный D packet: derive минимальную
complete final-audit fixture с field inventory и independent RHS, явно
зафиксировать manufactured applicability и ограниченный вывод об assembly,
а перенос на five-channel baseline оставить отдельным утверждением.
Никакой новый SCBA solve этим
решением не разрешён. Выключить каналы для зелёного full-state claim нельзя.

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
| Complete final U/G oracle готов | Exact seams/field inventory есть, coupled independent RHS нет | DEV05 blocked | данных недостаточно |
| Four absolute prefactors сертифицированы | Owning units/formulas и partial primary LO chain | Geometry/tolerance/primary gaps, DEV07 blocked | данных недостаточно |
| Current PV targeted runtime прошёл | GREEN2 fixture4533b105…,13+198 assertions; historical GREEN1 failure сохранён | FullCI/SCBA/science0 | сохранить |
| SCBA/science выполнены | Source/design, paper hashes, четыре bounded functional runs | Convergence/discretization/experiment не установлены | данных недостаточно |

Независимый CP-DESIGN принял четыре bounded functional scopes;
два blocked паспорта сохраняются. Reviewer проверил literal arithmetic,
shape/units, applicability, fixed RHS mutations и статусы. Автор документа
не закрывает собственный CP-CODE. Отдельный I medium выполняет только
accepted scope с envelope соответствующего DEV03/04/06/08; PV ratio
predicate получил CP-DESIGN/CP-CODE и bounded GREEN2 evidence. DEV05 и DEV07 не READY.
Полная six-oracle readiness этим документом не устанавливается.
