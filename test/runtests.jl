using Munkres
using Base.Test
using Combinatorics # Required for permutations() on 0.5


"""
Solve assignment problem via brute force exploration over all permutations
"""
function brute_force_optima(cost_matrix)
    (n,m) = size(cost_matrix)
    if n != m
        error("Non-square cost matrix, I haven't gotten around to doing that yet!")
    end
    initial_permuation = 1:n
    best_permuation = 1:n
    min_cost = typemax(eltype(cost_matrix))
    for p in permutations(initial_permuation)
        current_cost = zero(eltype(cost_matrix))
        for i = 1:n
            current_cost += cost_matrix[initial_permuation[i],p[i]]
        end
        if current_cost <= min_cost
            min_cost = current_cost
            best_permuation = p
        end
    end
    return best_permuation
end


# Small matrix designed to hit all steps
tst = [1 2 3;
       2 4 6;
       3 6 9]
@test munkres(tst) == brute_force_optima(tst)


for i=1:10
    tst = rand(8,8)
    @test munkres(tst) == brute_force_optima(tst)
end

for i=1:10
    # Exponentially distributed test for pathological floating point truncation
    tst = 10.^(50*rand(8,8))
    @test munkres(tst) == brute_force_optima(map(BigFloat, tst))
end

tst = -eye(10)
@test munkres(tst) == brute_force_optima(tst)

#test non-square behaviour on 1 x n and n x 1 arrays

for i=1:10
    tst = rand(1,10)
    @test munkres(tst)[1] == indmin(tst)
end

for i=1:10
    tst = rand(10,1)
    @test findfirst(munkres(tst).==1) == indmin(tst)
end

#test against solutions from Yi Cao's matlab code
#(see http://mathworks.com/matlabcentral/fileexchange/20328-munkres-assignment-algorithm)

tst = [ 678.0 24.0 601.0 34.0 747.0 455.0 871.0 714.0 64.0 547.0
        729.0 1.0 160.0 183.0 695.0 161.0 38.0 790.0 120.0 251.0
        983.0 472.0 960.0 63.0 994.0 761.0 375.0 287.0 959.0 698.0
        681.0 307.0 627.0 6.0 265.0 296.0 605.0 696.0 977.0 283.0
        685.0 203.0 811.0 665.0 422.0 430.0 202.0 317.0 879.0 934.0
        372.0 481.0 839.0 569.0 869.0 898.0 10.0 50.0 376.0 277.0
        395.0 634.0 85.0 156.0 219.0 590.0 924.0 194.0 434.0 294.0
        126.0 132.0 505.0 119.0 637.0 903.0 593.0 487.0 816.0 978.0
        48.0 665.0 766.0 707.0 834.0 291.0 853.0 813.0 504.0 716.0
        702.0 507.0 562.0 685.0 628.0 767.0 708.0 124.0 962.0 442.0
        827.0 71.0 946.0 672.0 117.0 47.0 592.0 876.0 614.0 115.0
        881.0 807.0 375.0 748.0 417.0 958.0 898.0 830.0 448.0 740.0
        207.0 681.0 506.0 159.0 294.0 547.0 109.0 876.0 802.0 231.0
        307.0 125.0 837.0 460.0 467.0 419.0 577.0 687.0 904.0 445.0
        757.0 190.0 26.0 601.0 49.0 71.0 814.0 476.0 268.0 998.0
        87.0 858.0 870.0 73.0 11.0 855.0 51.0 750.0 317.0 942.0
        747.0 687.0 424.0 630.0 175.0 446.0 371.0 14.0 902.0 406.0
        307.0 893.0 131.0 664.0 689.0 482.0 607.0 369.0 283.0 98.0
        442.0 59.0 979.0 463.0 161.0 240.0 328.0 949.0 409.0 299.0
        116.0 313.0 820.0 929.0 347.0 958.0 918.0 806.0 999.0 400.0]

p = [9, 2, 0, 4, 0, 7, 0, 0, 1, 0, 6, 0, 0, 0, 3, 5, 8, 10, 0, 0]

@test munkres(tst) == p

tst = [351.0 300.0 662.0 337.0 384.0 935.0 650.0 446.0 57.0 295.0 168.0 398.0 810.0 561.0 458.0 395.0 188.0 352.0 297.0 602.0
 221.0 662.0 123.0 448.0 255.0 124.0 724.0 935.0 422.0 340.0 335.0 945.0 368.0 423.0 347.0 809.0 823.0 769.0 879.0 886.0
 746.0 625.0 661.0 40.0 541.0 704.0 376.0 324.0 406.0 979.0 327.0 98.0 163.0 742.0 975.0 758.0 646.0 19.0 461.0 836.0
 791.0 962.0 62.0 123.0 813.0 394.0 364.0 297.0 266.0 192.0 863.0 298.0 569.0 462.0 698.0 876.0 317.0 227.0 400.0 869.0
 198.0 145.0 512.0 26.0 352.0 813.0 138.0 43.0 309.0 58.0 868.0 345.0 665.0 538.0 874.0 198.0 118.0 537.0 706.0 787.0
 370.0 975.0 516.0 944.0 310.0 175.0 11.0 895.0 932.0 818.0 264.0 754.0 727.0 938.0 844.0 451.0 11.0 828.0 974.0 654.0
 300.0 506.0 732.0 960.0 303.0 286.0 826.0 454.0 284.0 399.0 179.0 842.0 396.0 966.0 979.0 151.0 548.0 386.0 609.0 480.0
 479.0 209.0 945.0 947.0 493.0 906.0 97.0 488.0 204.0 817.0 681.0 201.0 583.0 740.0 218.0 237.0 662.0 417.0 814.0 672.0
 170.0 471.0 203.0 371.0 218.0 310.0 52.0 882.0 742.0 553.0 704.0 366.0 16.0 774.0 715.0 951.0 810.0 10.0 755.0 455.0
 536.0 155.0 664.0 400.0 807.0 385.0 837.0 439.0 804.0 885.0 371.0 929.0 118.0 210.0 553.0 104.0 367.0 767.0 255.0 50.0]

p = [9, 6, 18, 3, 4, 17, 16, 7, 13, 20]

@test munkres(tst) == p
