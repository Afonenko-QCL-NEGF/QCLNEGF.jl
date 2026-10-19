# DEV02 — паспорта независимых oracle O08-01…06

Статус: **CP-DESIGN PASS для bounded scopes O08-01/02/04/06**
(2026-10-19, независимый reviewer high). O08-03/O08-05 остаются
**design_blocked**. Полная готовность шести oracle и полное завершение
DEV02 acceptance **не установлены**.

Reviewed candidate SHA256:
`e3f577a55ca215f4e6c33fdc6883f4f38dfa4413b343c092389aac6e88832cca`.
Независимый CP-DESIGN receipt SHA256:
`80d0281721d77ba425a373824f49225c4713683958c2f646bba99cb1c7f4503e`.
Ни один статус ниже не означает исполненный test или scientific acceptance.
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
не переключались/не обновлялись. Все Core locators ниже относятся к `4af6e94`.

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
| O08-02 | design_accepted | Independent sine anchor + constant incompatible source; assessment routing явно synthetic |
| O08-03 | design_blocked | Не завершён complete manufactured final-state oracle с независимыми U/G/n/candidate/history |
| O08-04 | design_accepted | Nonzero signed public current с SI и literal two-state coherence/truncation |
| O08-05 | design_blocked | Четыре absolute channel passports ещё не имеют законченной primary→geometry→SI→tolerance цепочки; LO derivation — локальный кандидат |
| O08-06 | design_accepted | Literal discrete direct/FFT map и отдельный center-only continuum/tail oracle |

`design_accepted` означает только прошедший независимый CP-DESIGN
bounded functional-regression scope. DEV03/04/06/08 готовы к отдельной
реализации в соответствующих envelopes; runtime acceptance ещё не установлена.
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
healthy `n=ones(4)`, incompatible `n=2ones(4)`. Direct internal carrier
явно manufactured; это не изменение production grid/density. Тогда
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
Для assessment routing нужен отдельный явно synthetic carrier: pattern
`test/numerics/acceptance_metadata_registry.jl::assessed_fixture` можно
переиспользовать как **форму**, не запускать его seed/build как oracle.
У него manufacturing SCBA classification и pass values других gates —
явные routing inputs, не измеренные scientific values; не сохранять их как
реальное scientific report. В test assertion label это обозначить.
На synthetic healthy routing input `ConvergenceReport.passed==true`;
подстановка **только actual nonneutral Poisson diagnostics** даёт false и
messages для трёх соответствующих gates. Проверить каждый named failure,
чтобы отказ по unrelated field не маскировал dropped neutrality gate.
Если implementation не может выделить carrier без solve, ограничить runtime
direct Poisson slice и оставить assessment subcase blocked.

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

**Nonzero literal fixture.** NE=Nk=Nb=2; `wE=[1/4,3/4]`,
`wk=[1/2,3/2]`, gs=2. E0=.1 eV, L0=10 nm. В каждом (e,m)

```text
Q=[2 i;-i 1], B=[1 i;-i 2], F=I/2, D=I,
a=[1 2;3 4] indexed (e,m),
G< = i*a*Q, G> = -i*D, GR = -i*(a*Q+D)/2,
A=a*Q+D, Σ< = i*F, Σ> = -i*B, ΣR=-i*(B+F)/2.
```

Это current functional fixture, не Dyson/SCBA solution. Все correlation
matrices PSD; Green and family Keldysh identities заданы явно. Offdiagonal
entries влияют на literal trace: `tr(B*aQ-FD)=6a-1`, то есть
`S=[5 11;17 23]`. Independent scalar expansion даёт

```text
sum_m wk*S = [19,43]
jE = [19,43]/π * J0/E0_eV
sum_e (E0_eV*wE[e])*jE[e] = 37/π * J0
e*boundary_flux = 37/π * J0.
```

J0 вычислять в expected path из literals `e=1.602176634e-19 C`,
`hbar=1.054571817e-34 J s`, `E0_J=.1*e J`, `L0=1e-8 m`:
`J0=e*E0_J/(hbar*L0²)`; не читать expected из `scales.J0`.
Сначала отдельно проверить сам scale against this SI expression, затем
оба public API против retained expected values. Использовать physical
weights `E0_eV*wE`, не только dimensionless wE. Таким образом ошибка
общего multiplier/sign/2π у двух API не сокращается в одном agreement.

**Coordinate reconstruction — отдельная fixture.** Nz=3,Nb=2,
uniform Δx=1, constant mass mr>0; `H(g,g+1)=-t`,
`t=hbar²/(2*m0*mr*E0_J*L0²)` independently (m0=9.1093837139e-31 kg).
Literal orthonormal columns

```text
φ1=(1,1,0)/sqrt(2), φ2=(1,-1,2)/sqrt(6), Φ=[φ1 φ2],
G< = i*Q = [2i -1;1 i], NE=Nk=1, wE=wk=1, gs=2.
Re[(ΦG<Φ†)[2,1]]=-1/sqrt(3)
Re[(ΦG<Φ†)[3,2]]=+1/sqrt(3)
expected bonds = [ +2t/sqrt(3), -2t/sqrt(3) ] * J0/π.
```

χ=Φ for Δx=1. Diagonal-only G< gives both bonds 0; truncation to φ1
with retained G11=2i also gives 0. Это precisely derived effect данного
projected manufactured state; uniform J не требуется. Не нормировать
truncated density повторно и не сравнивать её с исходной как той же моделью.

**Норма/допуск/controls.** Componentwise absolute+relative comparison по
existing finite-chain current regression budget `atol=2e-14,rtol=2e-13`
применяется **после деления actual/expected current на независимый J0**,
а resolved — J0/E0_eV. Это existing small-functional arithmetic budget,
не 2e-14 A/m². Около нуля absolute branch применяется к scaled value;
weights перечислены выше. Expected scalar sums — dyadic exact до π/scale;
radicals и SI conversion ограничены этим regression budget, новый строгий
roundoff bound не заявлен. Mutations: sign, omitted factor2π, drop one wk,
drop wE on integration, diagonal-only or truncated reconstruction должны
нарушать соответствующий retained expected assertion.

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
| Keldysh/current/PV получают различающие RHS | Literal matrix/component sums и analytic PV omitted-cell/tail equations выше | CP-DESIGN PASS для bounded scopes, execution0 | сохранить |
| Poisson incompatibility видна через ζ/rP/rneutral | Constant RHS bordered equations, existing sine derivative anchor | Assessment carrier synthetic, full coupled state не доказан | сохранить |
| Complete final U/G oracle готов | Exact seams/field inventory есть, coupled independent RHS нет | DEV05 blocked | данных недостаточно |
| Four absolute prefactors сертифицированы | Owning units/formulas и partial primary LO chain | Geometry/tolerance/primary gaps, DEV07 blocked | данных недостаточно |
| Тесты/SCBA/science выполнены | Только source/design и paper hashes | Native/runtime/convergence/discretization/experiment not_measured | данных недостаточно |

Независимый CP-DESIGN принял четыре bounded functional scopes;
два blocked паспорта сохраняются. Reviewer проверил literal arithmetic,
shape/units, applicability, fixed RHS mutations и статусы. Автор документа
не закрывает собственный CP-CODE. Отдельный I medium выполняет только
accepted scope с envelope соответствующего DEV03/04/06/08; DEV05 и DEV07
не READY. Полная six-oracle readiness этим документом не устанавливается.
