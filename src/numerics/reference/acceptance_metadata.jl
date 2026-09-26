"""Version of the descriptive acceptance contract; formulas are evaluated elsewhere."""
const _ACCEPTANCE_METADATA_VERSION = "qcl-negf-acceptance-metadata-v1"

# A description is not a second implementation of a physical residual. Source
# references below identify the single evaluator; completeness is tested
# against every threshold emitted by _solution_validation_limits.
function _acceptance_metric_description(name::Symbol)
    label = String(name)
    validation = "src/numerics/reference/validation.jl"
    observable = "src/physics/observables.jl"
    audit = "src/numerics/reference/final_audit.jl"
    plain = "No quadrature: unweighted norms of the stored discrete arrays or block maxima."
    phase = "sum_e,m wE[e]*wk[m]; wE and wk are the declared energy and annular radial weights."
    none = (kind = "none", value = nothing, units = "dimensionless")
    additive = (
        kind = "additive",
        value = 1e-14,
        units = "scaled solver quantity in the denominator",
    )
    number = "N(G<)=gs/(2*pi)*sum_e,m wE[e]*wk[m]*real(-i*trace(G<[e,m]))"
    if name in (:r_D, :r_D_candidate)
        return (;
            formula = "max_e,m ||D*GR-I||F/(||D||F*||GR||F+sqrt(Nb)); D=E*I-H(U)-SigmaR",
            normalization = "Blockwise normwise backward residual; I and H use the declared orthonormal basis.",
            floor = (
                kind = "additive_identity_norm",
                value = "sqrt(Nb)",
                units = "dimensionless",
            ),
            weights = plain,
            source = validation*"::_dyson_residual",
            error_class = "algebraic",
            context = name===:r_D ?
                      "Stored total self-energy, published Hartree and Green state." :
                      "Fresh total candidate self-energy, same published Hartree and Green state.",
            tolerance = name===:r_D ? "r_D" : "r_K",
            applicability = "all_states",
        )
    elseif name in (:r_A, :r_A_candidate)
        return (;
            formula = "sqrt(sum_e,m ||A-GR*Gamma*GA||F^2)/(sqrt(sum_e,m ||A||F^2)+1e-14); Gamma=i*(SigmaGreater-SigmaLesser)",
            normalization = "Global unweighted Frobenius norm of the spectral array.",
            floor = additive,
            weights = plain,
            source = validation*"::_spectral_residual",
            error_class = "algebraic",
            context = name===:r_A ? "Stored total broadening." :
                      "Fresh total candidate broadening.",
            tolerance = name===:r_A ? "r_A" : "r_K",
            applicability = "all_states",
        )
    elseif name===:r_K
        return (;
            formula = "||G< - GR*SigmaCandidate<*GA||F/(||G<||F+1e-14)",
            normalization = "Global unweighted Frobenius norm; predicted G< is not divided by a fitted lambda.",
            floor = additive,
            weights = plain,
            source = validation*"::_keldysh_residual; "*audit*"::fresh_fixed_point_metrics",
            error_class = "nonlinear",
            context = "Fresh unmixed candidate, same accepted U and G.",
            tolerance = "r_K",
            applicability = "fresh_candidate",
        )
    elseif name===:r_Σ
        return (;
            formula = "max_family,component ||SigmaCandidate-SigmaStored||F/(||SigmaCandidate||F+1e-14)",
            normalization = "Max over retarded/lesser/greater and enabled scattering, total embedding, plus and minus embedding families.",
            floor = additive,
            weights = plain,
            source = validation*"::_selfenergy_residual; "*audit*"::fresh_fixed_point_metrics",
            error_class = "nonlinear",
            context = "Raw difference before mixing; no alpha factor.",
            tolerance = "r_Σ",
            applicability = "fresh_candidate",
        )
    elseif name in (:r_λ, :r_lambda_candidate)
        return (;
            formula = "abs(lambda-1); lambda=N(Graw<)/NDsheet; "*number,
            normalization = "NDsheet is the positive scaled target donor sheet density; raw G< is never normalized in this check.",
            floor = none,
            weights = phase*" Multiply by gs/(2*pi).",
            source = audit*"::fresh_fixed_point_metrics, charge_constraint_diagnostics",
            error_class = "nonlinear",
            context = name===:r_λ ? "Graw<=GR*SigmaStored<*GA." :
                      "Graw<=GR*SigmaCandidate<*GA.",
            tolerance = "r_λ",
            applicability = "fresh_candidate",
        )
    elseif name===:r_roundoff
        return (;
            formula = "max_enabled_channel ||Lambda1-Lambda2||F/(||Lambda1||F+1e-14)",
            normalization = "Norm of the real retarded reconstruction, per mechanism.",
            floor = additive,
            weights = "Hilbert quadrature uses wE; direct/product integration reverse the accumulation order; FFT repeats at twice valid padding.",
            source = "src/numerics/reference/kernels.jl::direct_hilbert_transform, product_integration_hilbert_transform; src/numerics/optimized/production_operators.jl::fft_hilbert_transform",
            error_class = "roundoff_indicator",
            context = "Empirical ordering/padding sensitivity, not a forward-error bound or finite-window error.",
            tolerance = "r_roundoff",
            applicability = "retarded_scattering_reconstruction",
        )
    elseif name in (:r_Jchange_scba, :r_Jchange)
        return (;
            formula = name===:r_Jchange_scba ?
                      "abs(Jnu-Jprev)/max(abs(Jnu),abs(Jprev),1e-12)" :
                      "max(last_outer.r_Jchange,abs(Jfinal-Jlast_outer)/max(abs(Jfinal),abs(Jlast_outer),1e-12))",
            normalization = "Relative change of the plus-boundary current; first comparison is unavailable.",
            floor = (kind = "maximum", value = 1e-12, units = "A/m^2"),
            weights = phase*" Current includes gs/(2*pi) and J0.",
            source = name===:r_Jchange_scba ?
                     "src/numerics/optimized/production_iteration.jl::_production_observable_measurement; src/numerics/reference/scba.jl::solve_scba" :
                     validation*"::validate_solution",
            error_class = "nonlinear",
            context = name===:r_Jchange_scba ? "Last accepted inner iteration record." :
                      "Outer history and independent final state comparison.",
            tolerance = "r_obs",
            applicability = name===:r_Jchange_scba ? "inner_history" : "outer_history",
        )
    elseif name in (:r_population_scba, :r_population)
        return (;
            formula = name===:r_population_scba ? "||pnu-pprev||2/(||pnu||2+1e-14)" :
                      "max(last_outer.r_population,||pfinal-plast_outer||2/(||pfinal||2+1e-14))",
            normalization = "p=real(diag(sheet_density_matrix))/NDsheet; vector norm, no renormalization of each level.",
            floor = (
                kind = "additive",
                value = 1e-14,
                units = "dimensionless population fraction",
            ),
            weights = phase*" Density matrix includes gs/(2*pi).",
            source = validation*"::validate_solution; src/numerics/optimized/production_iteration.jl::_production_observable_measurement",
            error_class = "nonlinear",
            context = "Observable stability is independent of fixed-point and conservation checks.",
            tolerance = "r_obs",
            applicability = name===:r_population_scba ? "inner_history" : "outer_history",
        )
    elseif name===:r_PSD
        return (;
            formula = "max_X,e,m max(0,-lambda_min(Herm(X))-b)/(max(||Herm(X)||2,floatmin(Float64))); b=64*Nb*eps(Float64)*max(||Herm(X)||2,1)",
            normalization = "X in {A,Gn_normalized,Gp_normalized,GammaTotal,Gn_raw,Gp_raw}; Herm(X)=(X+Xadjoint)/2 is used only to measure eigenvalues, never to repair the state.",
            floor = (
                kind = "maximum",
                value = floatmin(Float64),
                units = "scaled norm of X",
            ),
            weights = plain,
            source = validation*"::_negative_part, _green_positivity; src/numerics/reference/positivity_diagnostics.jl::psd_error_budget",
            error_class = "psd_backward_budget",
            context = "Final reference eigensolve; PSD_METRIC_VERSION=$(PSD_METRIC_VERSION); raw and number-constrained occupied/empty matrices are both checked; construction_scale=1 solver unit and construction_error=0.",
            tolerance = "r_PSD",
            applicability = "all_states",
        )
    elseif name in (:r_caus, :r_caus_sc, :r_caus_embedding)
        return (;
            formula = "max_e,m max(0,lambda_max((SigmaR-SigmaR_adjoint)/(2i)))/(||SigmaR||2+1e-14)",
            normalization = "Block spectral norm of the selected retarded self-energy.",
            floor = additive,
            weights = plain,
            source = validation*"::_positive_imaginary_part, _causality_residual",
            error_class = "algebraic",
            context = name===:r_caus ? "Total scattering plus embedding." :
                      name===:r_caus_sc ? "Internal scattering sum." : "Embedding sum.",
            tolerance = "r_caus",
            applicability = "all_states",
        )
    elseif name in (
        :r_selfenergy_identity,
        :r_selfenergy_identity_sc,
        :r_selfenergy_identity_embedding,
    )
        return (;
            formula = "max_e,m ||SigmaR-SigmaR_adjoint-(SigmaGreater-SigmaLesser)||F/(||SigmaR||F+||SigmaLesser||F+||SigmaGreater||F+1e-14)",
            normalization = "Sum of block Frobenius norms, not a spectral norm.",
            floor = additive,
            weights = plain,
            source = validation*"::_selfenergy_identity_residual",
            error_class = "algebraic",
            context = name===:r_selfenergy_identity ? "Total scattering plus embedding." :
                      name===:r_selfenergy_identity_sc ? "Internal scattering sum." :
                      "Embedding sum.",
            tolerance = "r_caus",
            applicability = "all_states",
        )
    elseif name in (
        :r_antihermitian_G,
        :r_antihermitian_Σ,
        :r_antihermitian_Σsc,
        :r_antihermitian_Σembedding,
    )
        return (;
            formula = "max_lesser,greater,e,m ||X+Xadjoint||F/(||X||F+1e-14)",
            normalization = "Frobenius norm of each lesser/greater block.",
            floor = additive,
            weights = plain,
            source = validation*"::_antihermiticity_residual",
            error_class = "algebraic",
            context = name===:r_antihermitian_G ? "Green components." :
                      name===:r_antihermitian_Σ ? "Total self-energy components." :
                      name===:r_antihermitian_Σsc ? "Scattering components." :
                      "Embedding components.",
            tolerance = "r_caus",
            applicability = "all_states",
        )
    elseif name===:r_sum
        return (;
            formula = "max_m ||sum_e wE[e]*A[e,m]/(2*pi)-I||F/sqrt(Nb)",
            normalization = "Frobenius matrix deficit divided by sqrt(Nb); not a missing-charge fraction or a current error.",
            floor = none,
            weights = "Energy weights wE/(2*pi), no radial weighting: worst represented k block.",
            source = validation*"::_spectral_sum_residual",
            error_class = "representation",
            context = "Finite energy window, orthonormal finite basis; full-axis identity is the comparison target.",
            tolerance = "r_sum",
            applicability = "all_states",
        )
    elseif name===:r_J
        return (;
            formula = "abs(real(PhiPlus+PhiMinus))/max(abs(real(PhiPlus)),abs(real(PhiMinus)),1e-12)",
            normalization = "Independent signed outgoing embedding number fluxes in scaled solver units.",
            floor = (
                kind = "maximum",
                value = 1e-12,
                units = "scaled boundary number flux",
            ),
            weights = phase*" Phi=gs/(2*pi)*sum weights*trace(SigmaGreater*G< - SigmaLesser*G>).",
            source = validation*"::validate_solution; "*observable*"::_boundary_flux_complex",
            error_class = "conservation",
            context = "Field-periodic embedding: plus/minus flux continuity; not the difference between successive iterates.",
            tolerance = "r_J",
            applicability = "field_periodic",
        )
    elseif name===:r_ImJ || startswith(label, "r_ImC_")
        return (;
            formula = "q=max_e,m abs(imag(trace_term)); s=max_e,m abs(real(trace_term)); s>1e-12 ? q/s : 100*q",
            normalization = "Pointwise imaginary cancellation before energy/momentum integration; current takes the maximum of plus/minus families.",
            floor = (
                kind = "branch_switch",
                value = 1e-12,
                units = "scaled pointwise trace; absolute branch multiplied by 100",
            ),
            weights = "No quadrature in this metric; pointwise maximum over all E,k.",
            source = validation*"::_pointwise_current_imaginary, _pointwise_collision_imaginary",
            error_class = "roundoff_indicator",
            context = name===:r_ImJ ?
                      "trace_term=Tr(Sigma>*G< - Sigma<*G>) for embedding." :
                      "trace_term=Tr(Sigma<*G> - Sigma>*G<) for named internal mechanism.",
            tolerance = "r_imag",
            applicability = name===:r_ImJ ? "field_periodic" : "named_scattering",
        )
    elseif name in (:r_roundJ, :r_roundP) || startswith(label, "r_roundC_")
        return (;
            formula = "abs(F)>1e-12 ? abs(F-Freverse)/abs(F) : 100*abs(F-Freverse)",
            normalization = "Maximum over plus/minus fluxes or internal channel powers; named collision uses its unprefactored signed integral.",
            floor = (
                kind = "branch_switch",
                value = 1e-12,
                units = "scaled integrated flux/power; absolute branch multiplied by 100",
            ),
            weights = phase*" Power additionally uses E[e]. Forward order sorts by (abs(E),E); reverse inverts that order.",
            source = validation*"::_roundoff_metric; "*observable*"::collision_balance, power_balance",
            error_class = "roundoff_indicator",
            context = "Change of accumulation order is empirical roundoff sensitivity, not an error bound.",
            tolerance = "r_roundoff",
            applicability = name===:r_roundJ ? "field_periodic" :
                            name===:r_roundP ? "internal_scattering" : "named_scattering",
        )
    elseif startswith(label, "r_C_")
        return (;
            formula = "c*abs(real(sum_e,m wE*wk*(incoming-outgoing)))/(c*sum_e,m wE*wk*(abs(incoming)+abs(outgoing))+1e-14); c=gs/(2*pi)",
            normalization = "incoming=Tr(Sigma<*G>), outgoing=Tr(Sigma>*G<); signed integral before absolute value.",
            floor = (
                kind = "additive",
                value = 1e-14,
                units = "scaled collision number rate",
            ),
            weights = phase,
            source = observable*"::_collision_sum, collision_balance",
            error_class = "conservation",
            context = "Electron-number conservation for the named internal channel; phonon energy exchange need not vanish away from equilibrium.",
            tolerance = "r_C",
            applicability = "named_scattering",
        )
    elseif name===:r_power
        return (;
            formula = "abs(sum_channels Pchannel+J*Vperiod)/(sum_channels abs(Pchannel)+abs(J*Vperiod)+1e-14*P0)",
            normalization = "Pchannel is electron energy gain in W/m^2; J is charge times plus-boundary number flux, Vperiod=Fbias*Lp.",
            floor = (kind = "additive", value = "1e-14*P0", units = "W/m^2"),
            weights = phase*" Collision powers include gs/(2*pi), E[e] and P0.",
            source = observable*"::_collision_power_bar, power_balance",
            error_class = "conservation",
            context = "Stationary field-periodic electron system with declared internal channels. An open finite contact problem requires its contact energy flux terms; this formula is not that model.",
            tolerance = "r_power",
            applicability = "field_periodic",
        )
    elseif name===:r_ImP
        return (;
            formula = "100*max_internal_channel abs(imag(Pchannel_bar))",
            normalization = "Absolute imaginary part of dimensionless internal energy exchange, multiplied by 100; not a relative power error.",
            floor = none,
            weights = phase*" Multiply integrand by E[e] and gs/(2*pi).",
            source = observable*"::power_balance; "*validation*"::validate_solution",
            error_class = "roundoff_indicator",
            context = "Pchannel_bar=Pchannel/P0 before discarding its imaginary part.",
            tolerance = "r_imag",
            applicability = "internal_scattering",
        )
    elseif name===:r_P
        return (;
            formula = "||L*Ucandidate-b||2/(||b||2+lambda_P*1e-14); b=-lambda_P*(ND-n)",
            normalization = "Discrete periodic Poisson residual at the independently recomputed gauge-fixed potential.",
            floor = (
                kind = "additive",
                value = "lambda_P*1e-14",
                units = "scaled Poisson source",
            ),
            weights = "Euclidean algebraic vector norm; L uses cell spacing wZ[1] and harmonic interface permittivities.",
            source = "src/numerics/reference/poisson.jl::solve_periodic_poisson",
            error_class = "algebraic",
            context = "Neutral field-periodic source and bordered mean-zero gauge.",
            tolerance = "r_P",
            applicability = "field_periodic",
        )
    elseif name===:r_U
        return (;
            formula = "max_z abs(Ucandidate[z]-Upublished[z])",
            normalization = "Absolute Hartree energy change expressed in E0, before Poisson mixing; physical energy error indicator = r_U*E0.",
            floor = none,
            weights = "Spatial maximum, no quadrature.",
            source = validation*"::validate_solution; src/numerics/reference/poisson.jl::solve_periodic_poisson",
            error_class = "nonlinear",
            context = "Independent Poisson map of the same published density.",
            tolerance = "r_U",
            applicability = "field_periodic",
        )
    elseif name===:r_n
        return (;
            formula = "max(last_outer.r_n,sqrt(sum_z wZ[z]*(nfinal[z]-nlast_outer[z])^2)/(NDsheet/sqrt(Lp)))",
            normalization = "Weighted L2 density change divided by NDsheet/sqrt(Lp); last_outer.r_n uses the preceding outer density by the same formula.",
            floor = none,
            weights = "Spatial quadrature wZ; no fitted density renormalizer in the residual.",
            source = validation*"::validate_solution; src/numerics/reference/solver.jl::solve",
            error_class = "nonlinear",
            context = "Scaled density n*L0^3, Lp/L0 and sheet density ND*L0^2; requires an outer history.",
            tolerance = "r_n",
            applicability = "outer_history",
        )
    elseif name===:r_neutral
        return (;
            formula = "abs(sum_z wZ[z]*(ND[z]-n[z]))/sum_z wZ[z]*ND[z]",
            normalization = "Positive target donor sheet density, with no additive floor.",
            floor = none,
            weights = "Spatial quadrature wZ.",
            source = "src/numerics/reference/poisson.jl::solve_periodic_poisson",
            error_class = "conservation",
            context = "Neutrality of the stored density; normalization alone is not a fresh Keldysh certificate.",
            tolerance = "r_neutral",
            applicability = "field_periodic",
        )
    elseif name===:r_ζ
        return (;
            formula = "abs(zeta), from [L/sP v; v^T 0]*[U;zeta]=[b/sP;0]; v=ones(Nz)/sqrt(Nz), sP=max_i sum_j abs(Lij)",
            normalization = "Absolute Lagrange multiplier of the scaled bordered Poisson system; not an independently fitted electrostatic offset.",
            floor = none,
            weights = "Unweighted mean-zero gauge; matrix row-sum scaling sP.",
            source = "src/numerics/reference/poisson.jl::solve_periodic_poisson",
            error_class = "algebraic",
            context = "Periodic Poisson null mode and source compatibility.",
            tolerance = "r_ζ",
            applicability = "field_periodic",
        )
    elseif name in (:r_tail_low, :r_tail_high, :r_tail_k)
        return (;
            formula = "abs(sum_selected p[e,m])/sum_all p[e,m]; p=wE[e]*wk[m]*real(-i*trace(G<[e,m]))",
            normalization = "Total positive occupied discrete weight; nonfinite or nonpositive total returns Inf, never a floored denominator.",
            floor = none,
            weights = phase*" The common gs/(2*pi) factor cancels.",
            source = validation*"::_edge_metrics",
            error_class = "representation",
            context = name===:r_tail_k ?
                      "Last max(1,ceil(momentum_tail_window_fraction*Nk)) radial nodes." :
                      "First/last max(1,ceil(energy_tail_window_fraction*NE)) energy nodes, as selected by the metric suffix.",
            tolerance = "r_tail",
            applicability = "all_states",
        )
    elseif name===:r_edge_gamma
        return (;
            formula = "max_Eendpoint,m ||GammaSc[E,m]||F / max_all_E,m ||GammaSc[E,m]||F",
            normalization = "Internal scattering broadening only; zero denominator gives zero in the evaluator and is inapplicable when no internal channels exist.",
            floor = none,
            weights = "Unweighted endpoint/block maxima; no embedding broadening included.",
            source = validation*"::_edge_metrics",
            error_class = "representation",
            context = "Finite-window convergence heuristic, not a universal identity for all baths or boundary conditions.",
            tolerance = "r_edge",
            applicability = "internal_scattering",
        )
    elseif name===:r_edge_spectral
        return (;
            formula = "max_Eendpoint,m abs(real(trace(A[E,m])))/max_all_E,m abs(real(trace(A[E,m])))",
            normalization = "Endpoint spectral trace relative to its global maximum; zero maximum returns Inf.",
            floor = none,
            weights = "Unweighted endpoint/node maxima; no energy integration.",
            source = validation*"::_edge_metrics",
            error_class = "representation",
            context = "Endpoint heuristic assessed together with integrated spectral sum and independent energy-window/grid tests.",
            tolerance = "r_edge",
            applicability = "all_states",
        )
    elseif name===:r_energy_electron_electron
        return (;
            formula = "abs(sum_e,m wE*wk*E*(real(incoming)-real(outgoing)))/max(sum_e,m wE*wk*abs(E)*(abs(real(incoming))+abs(real(outgoing))),floatmin(Float64))",
            normalization = "Electronic energy moment of the declared e-e/SPPA family; real parts are taken before accumulation.",
            floor = (
                kind = "maximum",
                value = floatmin(Float64),
                units = "scaled electronic collision energy",
            ),
            weights = phase*" Include signed E in numerator and abs(E) in activity denominator.",
            source = "src/numerics/optimized/electron_electron.jl::electron_electron_collision_diagnostics",
            error_class = "conservation",
            context = "Zero electronic energy transfer applies only to a conserving isolated electron interaction. The prescribed SPPA bath reports this diagnostic without a zero-transfer gate; its power remains in the total power balance.",
            tolerance = "r_power",
            applicability = "electron_electron_energy_conserving_closure",
        )
    elseif name===:r_hot_lo_power
        return (;
            formula = "abs(Pelectron_collision+Pelectron_to_LO_events)/max(abs(Pelectron_collision),abs(Pelectron_to_LO_events),1e-30)",
            normalization = "Independent LO event-count energy transfer versus published electron collision energy, both W/m^2.",
            floor = (kind = "maximum", value = 1e-30, units = "W/m^2"),
            weights = phase*" Collision moment includes E; event count uses the fixed LO quantum.",
            source = audit*"::_lo_balance_diagnostics",
            error_class = "conservation",
            context = "LO enabled with lo_population=rate_balance; thermal fixed bath does not use this extra gate.",
            tolerance = "r_power",
            applicability = "hot_lo",
        )
    end
    throw(ArgumentError("acceptance metadata is undefined for metric $name"))
end

function _acceptance_metric_metadata(name::Symbol, solution::NEGFSolution)
    description=_acceptance_metric_description(name)
    problem=solution.problem
    active=problem.kernels.enabled
    scope=description.applicability
    applicable, reason=true, nothing
    if scope=="named_scattering"
        label=String(name)
        prefix=startswith(label, "r_roundC_") ? "r_roundC_" :
               startswith(label, "r_ImC_") ? "r_ImC_" : "r_C_"
        mechanism=Symbol(chopprefix(label, prefix))
        applicable=mechanism in active
        reason=applicable ?
               "Named electron-number-conserving internal channel is enabled." :
               "Named scattering channel is disabled."
    elseif scope in ("internal_scattering", "retarded_scattering_reconstruction")
        applicable=!isempty(active)
        reason=applicable ? "At least one internal scattering channel is enabled." :
               "No internal scattering channels are enabled."
        algorithms=get(
            get(solution.observables, :restart_contract, Dict()),
            "algorithms",
            Dict(),
        )
        if scope=="retarded_scattering_reconstruction" &&
           get(algorithms, "retarded_real_part", nothing)=="drop"
            applicable=false
            reason="The declared approximation drops the real retarded reconstruction."
        end
    elseif scope=="electron_electron"
        applicable=:electron_electron in active
        reason=applicable ?
               "The declared e-e/SPPA channel is enabled; the effective-bath qualification applies." :
               "The e-e/SPPA channel is disabled."
    elseif scope=="hot_lo"
        applicable=problem.scattering.LO && problem.models.lo_population===:rate_balance
        reason=applicable ? "LO rate-balance population is enabled." :
               "The LO rate-balance population model is not enabled."
    else
        reason=scope=="outer_history" ?
               "Coupled Poisson-SCBA acceptance requires outer history; absent history is not a pass." :
               scope=="inner_history" ?
               "An inner iteration sequence is required; absent history is not a pass." :
               scope=="fresh_candidate" ?
               "Fresh unmixed candidate is required; unavailable audit is not a pass." :
               scope=="field_periodic" ?
               "Current implementation is the stationary field-periodic embedding model." :
               "Applies to every represented state in the declared orthonormal basis."
    end
    metadata=Dict{String,Any}(
        "registry_version"=>_ACCEPTANCE_METADATA_VERSION,
        "metric_id"=>String(name),
        "theory_reference"=>"docs/src/theory/12_validation.md#eq-validation-residuals; docs/src/theory/11_observables.md",
        "formula"=>description.formula,
        "formula_source"=>description.source,
        "formula_context"=>description.context,
        "units"=>"dimensionless",
        "normalization"=>description.normalization,
        "denominator_floor"=>Dict(String(k)=>v for (k, v) in pairs(description.floor)),
        "weights"=>description.weights,
        "applicability"=>Dict("applies"=>applicable, "scope"=>scope, "reason"=>reason),
        "threshold_origin"=>Dict(
            "kind"=>"project_acceptance_policy",
            "parameter"=>"SolverTolerances."*description.tolerance,
            "source"=>"src/domain/types.jl::SolverTolerances; frozen resolved configuration",
            "literature_threshold"=>false,
            "interpretation"=>"Acceptance target, not an estimated error or a universal NEGF constant.",
        ),
        "error_budget"=>Dict(
            "class"=>description.error_class,
            "estimated_bound"=>nothing,
            "bound_status"=>"not_measured",
            "interpretation"=>description.error_class=="psd_backward_budget" ?
                              "Declared 64*Nb*eps construction/eigensolve allowance is applied before residual; it is not a spectral discretization budget." :
                              description.error_class=="roundoff_indicator" ?
                              "Empirical arithmetic sensitivity; no rigorous forward-error bound is inferred." :
                              description.error_class=="nonlinear" ?
                              "Raw residual; observable error requires conditioning/contraction evidence and a strict final audit." :
                              description.error_class=="representation" ?
                              "Domain/resolution indicator; window, step, k cutoff and basis errors must be measured independently." :
                              "Residual of the declared discrete model; it does not bound the observable error without an independent refinement study.",
            "discretization_evidence"=>"campaign-level grid comparisons; not certified by this single-state metric",
        ),
    )
    algorithms=get(
        get(solution.observables, :restart_contract, Dict()),
        "algorithms",
        Dict(),
    )
    if name===:r_roundoff && get(algorithms, "hilbert", nothing)=="product_integration"
        metadata["formula"]="max_enabled_channel ||Lambda1-Lambda2||F/max(||Lambda1||F,eps(Float64))"
        metadata["denominator_floor"]=Dict(
            "kind"=>"maximum",
            "value"=>eps(Float64),
            "units"=>"scaled real retarded self-energy norm",
        )
        metadata["formula_context"]="Product integration compares opposite matrix-column orders; its denominator differs from direct/FFT reconstruction."
    end
    return metadata
end
