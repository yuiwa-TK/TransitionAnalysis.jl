using Test
using TransitionAnalysis

@testset "TransitionAnalysis.jl" begin
    # Write your tests here.
    @testset "Constants" begin
        cc = TransitionAnalysis.Consts
        @test cc.γ==1.4
        @test cc.pr==0.72
    end;

    @testset "Basefunc" begin
        @test TransitionAnalysis.Sutherland(1.0)==1.0
    end;

    @testset "EmpericalRelation" begin
        cc = TransitionAnalysis.Consts
        r  = cc.pr^(1/3)
        γ  = cc.γ
        TAemp = TransitionAnalysis.EmpericalRelation
        @test TAemp.cf_lami__incompressible(100)==0.0664
        @test TAemp.recovery_temp_laminar(0.1)== 1.0 + 0.5*r*(γ-1)*0.1*0.1
    end;
end;
