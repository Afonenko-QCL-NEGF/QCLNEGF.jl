module PublicBasisDiagnosticsTests
include("../support/common.jl")
include("../support/public_state_fixture.jl")

@testset "Public basis subset keeps conditional pzp thresholds and does not build a transport state" begin
    available=isdefined(QCLNEGF,:basis_validation_diagnostics)
    @test available
    if available
        f=public_state_fixture()
        basis=f.problem.basis
        changes=(r_orth=0.0,r_eigen=0.0,r_translation=0.02,r_translation_nearest=0.02,
            r_translation_two_pairs=0.01,T₊=zeros(ComplexF64,1,1),T₋=zeros(ComplexF64,1,1),localization=:pzp)
        controlled=typeof(basis)((get(changes,name,getfield(basis,name)) for name in fieldnames(typeof(basis)))...)
        h=zeros(ComplexF64,f.problem.numerical.N_k,1,1)
        small=basis_validation_diagnostics(controlled,h,f.problem.numerical)
        @test small.report.passed
        @test small.report.metrics[:translation_two_pairs]==0.01
        @test small.limits[:translation_two_pairs]==0.1
        n=NumericalParameters(; (name=>(name===:P_basis ? 3 : getfield(f.problem.numerical,name))
            for name in fieldnames(NumericalParameters))...)
        larger=basis_validation_diagnostics(controlled,h,n)
        @test !larger.report.passed && larger.limits[:translation_two_pairs]==0.005
        plain=typeof(basis)((name===:localization ? :none : getfield(controlled,name)
            for name in fieldnames(typeof(basis)))...)
        @test !haskey(basis_validation_diagnostics(plain,h,n).report.metrics,:translation)
        @test_throws DimensionMismatch basis_validation_diagnostics(controlled,zeros(ComplexF64,1,2,2),n)
        small.limits[:basis_eigen]=9.0
        @test basis_validation_diagnostics(controlled,h,f.problem.numerical).limits[:basis_eigen]==1e-12
    end
end

@testset "Projected spectrum enforces Hermiticity before an independent two-level eigen oracle" begin
    available=isdefined(QCLNEGF,:projected_hamiltonian_spectrum)
    @test available
    if available
        h=ComplexF64[2 im;-im 2]
        result=projected_hamiltonian_spectrum(h)
        @test result.eigenvalues≈[1.0,3.0] rtol=3e-15
        @test h*result.eigenvectors≈result.eigenvectors*Diagonal([1.0,3.0]) rtol=3e-15
        @test result.hermiticity_residual==0 && result.eigen_residual<1e-15
        result.eigenvectors[1,1]=99
        @test h==ComplexF64[2 im;-im 2]
        zeroed=projected_hamiltonian_spectrum(zeros(ComplexF64,2,2))
        @test zeroed.eigenvalues==[0.0,0.0] && zeroed.eigen_residual==0
        @test_throws DomainError projected_hamiltonian_spectrum(ComplexF64[1 im;im 2])
        @test_throws DomainError projected_hamiltonian_spectrum(fill(ComplexF64(NaN),2,2))
        @test_throws DimensionMismatch projected_hamiltonian_spectrum(zeros(ComplexF64,2,3))
        # Independent constant-mass finite-difference sine spectrum. Expected
        # eigenvalues do not call eigen or either projected-spectrum evaluator.
        count=12; step=0.2; coefficient=0.5; mass=0.1
        fd=build_bdd_hamiltonian(zeros(count),fill(mass,count),step,coefficient)
        expected=[2coefficient/(mass*step^2)*(1-cos(j*π/(count+1))) for j in 1:count]
        @test projected_hamiltonian_spectrum(fd).eigenvalues≈expected rtol=1e-12
    end
end
end
