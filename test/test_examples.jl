@eval module $(gensym())
using ITensorIMPSTools: ITensorIMPSTools
using Test: @test, @testset

@testset "Test examples" begin
  example_files = ["ising.jl"]
  @testset "Test $example_file" for example_file in example_files
    @test include(joinpath(pkgdir(ITensorIMPSTools), "examples", example_file)) == nothing
  end
end
end
