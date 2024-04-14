LoadPackage("QDistRnd");

F := GF(2);

CyclicShiftMat := function(l, m, identityFirst)
    local ml, mm;
    if identityFirst then
        ml := IdentityMat(l, F);
        mm := PermutationMat(CycleFromList([1..m]), m, F);
    else
        ml := PermutationMat(CycleFromList([1..l]), l, F);
        mm := IdentityMat(m, F);
    fi;
    return KroneckerProduct(ml, mm);
end;

# Determine if the given code has toric layout(Lemma4)
# If the code has toric layout, it must be single connected(Lemma3)
IsToricLayout := function(aGenerators, bGenerators, sizeM)
    local i, j;
    for i in aGenerators do
        for j in bGenerators do
            if (i[2] * j[2] = sizeM) and (Size(GroupByGenerators([i[1], j[1]])) = sizeM) then
                return true;
            fi;
        od;
    od;
    return false;
end;


GetSearchPool := function(l, m, gx, gy, x, y)
    local combs, selected, allMats, pool, i, j, elements, inverseElements, allPotentialGenerators, generators, generator;
    combs := Combinations(Concatenation(List([1..l], ii -> [0, ii-1]), List([1..m], ii -> [1, ii-1])), 3);
    
    allPotentialGenerators := [];
    for selected in combs do
        elements := List(selected, ii -> [gx, gy][ii[1] + 1] ^ ii[2]);
        inverseElements := List(elements, Inverse);
        generators := [];
        for i in [1..2] do
            for j in [i+1..3] do
                generator := elements[i] * inverseElements[j];
                Add(generators, [generator, Order(generator)]);
            od;
        od;
        Add(allPotentialGenerators, generators);
    od;

    pool := [];
    for i in [1..Length(combs)] do
        for j in [i..Length(combs)] do
            Add(pool, [i, j]);
        od;
    od;
    return MakeImmutable([combs, allPotentialGenerators, pool]);
end;


# Construct the parity check matrix from the given A and B blocks
ConstructParityCheckMatrix := function(A, B)
    local AT, BT, Hx, Hz;
    AT := TransposedMat(A);
    BT := TransposedMat(B);
    # Hx = [A | B]
    Hx := TransposedMat(Concatenation(AT, BT));
    # Hz = [BT | AT]
    Hz := TransposedMat(Concatenation(B, A));
    return [Hx, Hz];
end;

# SplitListIntoChunks := function(list, chunk_size)
#     local result, i, e;
#     result := [];
#     i := 0;
#     while i * chunk_size < Length(list) do
#         Add(result, list{[i * chunk_size + 1..Minimum((i + 1) * chunk_size, Length(list))]});
#         i := i + 1;
#     od;
#     return MakeReadOnlyObj(result);
# end;

ComputeK := function(x, y, combs, code)
    local A, B;
    A := Sum(List(combs[code[1]], ii -> [x, y][ii[1] + 1] ^ ii[2]));
    B := Sum(List(combs[code[2]], ii -> [x, y][ii[1] + 1] ^ ii[2]));
    return 2 * (Length(A) - Rank(Concatenation(A, B)));
end;

SatisfyEncodingRateBound := function(k, n, encodingRateLowerBound, targetK)
    if targetK > 0 and k <> targetK then
        return false;
    fi;
    return k / n / 2 >= encodingRateLowerBound;
end;

ComputeDistance := function(x, y, combs, code, numInformationSets, distanceLowerBound)
    local Hs, A, B;
    A := Sum(List(combs[code[1]], ii -> [x, y][ii[1] + 1] ^ ii[2]));
    B := Sum(List(combs[code[2]], ii -> [x, y][ii[1] + 1] ^ ii[2]));
    Hs := ConstructParityCheckMatrix(A, B);
    return DistRandCSS(Hs[1], Hs[2], numInformationSets, distanceLowerBound, 2: field:=F);
end;


# Given the parameter l and m, search the bivariate bicyclic code
SearchBBCodes := function(l, m, resultsFilepath, encodingRateLowerBound, distanceLowerBound, numInformationSets, targetK)
    local x_shift_mat, y_shift_mat, gs, ps, toric_layout_codes, codes_for_distance_calc, combs, tasks, allPotentialGenerators, satisfiedCodes, distances, d, code, i, ks, n, k;

    x_shift_mat := CyclicShiftMat(l, m, false);
    y_shift_mat := CyclicShiftMat(l, m, true);
    gs := GeneratorsOfGroup(AbelianGroup([l, m]));

    ps := GetSearchPool(l, m, gs[1], gs[2], x_shift_mat, y_shift_mat);
    combs := ps[1];
    # allMats := ps[2];
    allPotentialGenerators := ps[2];
    Print("Search pool size: ", Length(ps[3]), "\n");

    # Filter toric layout codes
    toric_layout_codes := Filtered(ps[3], ij -> IsToricLayout(allPotentialGenerators[ij[1]], allPotentialGenerators[ij[2]], l * m));
    Print("Number of codes with toric layout: ", Length(toric_layout_codes), "\n");

    # Filter encoding rate
    tasks := List(toric_layout_codes, code -> RunTask(ComputeK, x_shift_mat, y_shift_mat, combs, code));
    ks := List(tasks, TaskResult);
    codes_for_distance_calc := Filtered([1..Length(toric_layout_codes)], i -> SatisfyEncodingRateBound(ks[i], 2 * l * m, Rat(encodingRateLowerBound), targetK));
    Print("Number of code satisfying encoding rate conditions: ", Length(codes_for_distance_calc), "\n");

    # Compute the distance
    tasks := List(codes_for_distance_calc, i -> RunTask(ComputeDistance, x_shift_mat, y_shift_mat, combs, toric_layout_codes[i], numInformationSets, distanceLowerBound));
    distances := List(tasks, TaskResult);

    satisfiedCodes := [];
    n := 2 * l * m;
    for i in [1..Length(codes_for_distance_calc)] do
        d := distances[i];
        if d > 0 then
            code := toric_layout_codes[codes_for_distance_calc[i]];
            k := ks[codes_for_distance_calc[i]];
            Add(satisfiedCodes, rec(parameter:=[n, k, d], encodingRate:=k/n, As:=combs[code[1]], Bs:=combs[code[2]]));
        fi;
    od;
    Print("Number of code satisfying all conditions: ", Length(satisfiedCodes), "\n");

    if Length(satisfiedCodes) > 0 then
        PrintCSV(resultsFilepath, satisfiedCodes);
    fi;
end;


# numInformationSets := 10;
# chunkSize := 1000;
# resultsFilepath := "out/codes.csv";

# l := 12;
# m := 12;
# distanceLowerBound := 10;
# encodingRateLowerBound := 1 / 25;
# targetK := 0;
# SearchBBCodes(l, m, resultsFilepath, encodingRateLowerBound, distanceLowerBound, 10, 5000, targetK);

# l := 6;
# m := 6;
# distanceLowerBound := 5;
# encodingRateLowerBound := 1 / 13;
# targetK := 0;


# SearchBBCodes(l, m, resultsFilepath, encodingRateLowerBound, distanceLowerBound, numInformationSets, chunkSize, targetK);