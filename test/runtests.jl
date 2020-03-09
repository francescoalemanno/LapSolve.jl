using SolveLAP
using Test

@testset "StiffWrapper" begin
     A = [3.0 5.0 6.0 9.0 9.0;
          2.0 4.0 3.0 5.0 6.0;
          4.0 1.0 4.0 9.0 9.0;
          4.0 7.0 8.0 5.0 10.0;
          1.0 1.0 10.0 10.0 5.0]
     M = [2.0 4.0 5.0 8.0 8.0 18.0 Inf Inf Inf Inf;
          1.0 3.0 2.0 4.0 5.0 Inf 18.0 Inf Inf Inf;
          3.0 0.0 3.0 8.0 8.0 Inf Inf 18.0 Inf Inf;
          3.0 6.0 7.0 4.0 9.0 Inf Inf Inf 18.0 Inf;
          0.0 0.0 9.0 9.0 4.0 Inf Inf Inf Inf 18.0;
          18.0 Inf Inf Inf Inf 0.0 0.0 0.0 0.0 0.0;
          Inf 18.0 Inf Inf Inf 0.0 0.0 0.0 0.0 0.0;
          Inf Inf 18.0 Inf Inf 0.0 0.0 0.0 0.0 0.0;
          Inf Inf Inf 18.0 Inf 0.0 0.0 0.0 0.0 0.0;
          Inf Inf Inf Inf 18.0 0.0 0.0 0.0 0.0 0.0]
      W=SolveLAP.StiffWrapper(A,2)
      okay=true
      for i in eachindex(M,W)
            M[i]==W[i] || (okay=false; break)
      end
      @test typeof(first(eachindex(M,W))) <: CartesianIndex
      @test okay
end

@testset "simple examples" begin
      @testset "simple examples" begin
          A = [ 0.891171  0.0320582   0.564188  0.8999    0.620615;
                0.166402  0.861136    0.201398  0.911772  0.0796335;
                0.77272   0.782759    0.905982  0.800239  0.297333;
                0.561423  0.170607    0.615941  0.960503  0.981906;
                0.748248  0.00799335  0.554215  0.745299  0.42637]

          assign, cost = @inferred solve_lap(A)
          @test assign == [2, 3, 5, 1, 4]

          B = [ 24     1     8;
                 5     7    14;
                 6    13    20;
                12    19    21;
                18    25     2]

          assign, cost = solve_lap(B)
          @test assign == [2, 1, 0, 0, 3]
          @test cost == 8

          assign, cost = solve_lap(B')
          @test assign == [2, 1, 5]
          @test cost == 8

          assign, cost = solve_lap([ifelse(i-j==0,0,i*j) for i in 1:5, j in 1:5])
          @test assign == [1, 2, 3, 4, 5]
          @test cost == 0
      end
end



@testset "UInt16" begin
    M=UInt16[28092 44837 19882 39481 59139;
             26039 46258 38932 51057 9;
             11527 59487 61993 29072 8734;
             10691 16977 12796 16370 14266;
             5199  42319 34194 41332 16472]
    assign,cost=solve_lap(M)
    @test assign == [3, 5, 4, 2, 1]
    @test cost   == 71139
end

@testset "UInt8" begin
    M=UInt8[67  228 135 197 244;
            112 44  84  206 31;
            225 103 231 225 227;
            170 37  135 9   130;
            110 22  133 77  96]
    assign,cost=solve_lap(M)
    @test assign == [1, 5, 2, 4, 3]
    @test cost   == 343
end

@testset "stiff Problem" begin
    M=[0.08886586609282232 0.5001403306848344;
       0.6848622386100172 0.1379247179005112;
       0.9410463734407415 0.30774601623574327]
    A=[Inf Inf 1000.0;
     0.5266269523801999 Inf Inf;
      0.9238847081919377 Inf 0.6017336762941108]
    @test solve_stiff_lap(M)[1:3]==[(1, 1), (2, 2), (3, -1)]
    @test solve_stiff_lap(A)[1:3]== [(1, -1), (2, 1), (3, 3)]
    @test_throws ErrorException solve_lap(A)
end
