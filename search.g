LoadPackage("QDistRnd");

# set up the code in binary field
F := GF(2);


# Create the cyclic shift matrix
CyclicPermutationMat := function(s)
    local perm;
    perm := PermList(Concatenation([2..s], [1]));
    return PermutationMat(perm, s, F);
end;


# Create x/y cyclic matrix
CyclicMat := function(l, m, identityFirst)
    local ml, mm;
    if identityFirst then
        ml := IdentityMat(l, F);
        mm := CyclicPermutationMat(m);
    else
        ml := CyclicPermutationMat(l);
        mm := IdentityMat(m, F);
    fi;
    return KroneckerProduct(ml, mm);
end;


# Get direct product of CyclicGroup(l) and CyclicGroup(m)
# and its two generators
ConstructGroup := function(l, m)
    local cl, cm, M, gx, gy;
    cl := CyclicGroup(l);
    cm := CyclicGroup(m);
    M := DirectProduct(cl, cm);
    gx := M.1;
    gy := GeneratorsOfGroup(M)[Length(GeneratorsOfGroup(cl)) + 1];
    return [Size(M), gx, gy];
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


# Get the search pool for the bivariate bicyclic code
GetSearchPool := function(l, m, gx, gy, x, y)
    local combs, selected, allMats, pool, i, j, elements, inverseElements, allPotentialGenerators, generators, generator;
    combs := Combinations(Concatenation(List([1..l], ii -> [0, ii-1]), List([1..m], ii -> [1, ii-1])), 3);
    
    allMats := List(combs, selected -> Sum(List(selected, ii -> [x, y][ii[1] + 1] ^ ii[2])));

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
    return [combs, allMats, allPotentialGenerators, pool];
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

SplitListIntoChunks := function(list, chunk_size)
    local result, i, e;
    result := [];
    i := 0;
    while i * chunk_size < Length(list) do
        Add(result, list{[i * chunk_size + 1..Minimum((i + 1) * chunk_size, Length(list))]});
        i := i + 1;
    od;
    return result;
end;

CalculateFor := function(ijList, combs, allMats, allPotentialGenerators, sizeM, n, encodingRateLowerBound, numInformationSets, distanceLowerBound, targetK)
    local codes, ij, i, j, A, B, encodingRate, Hs, k, d;
    codes := [];
    for ij in ijList do
        i := ij[1];
        j := ij[2];
        if not IsToricLayout(allPotentialGenerators[i], allPotentialGenerators[j], sizeM) then
            continue;
        fi;
        A := allMats[i];
        B := allMats[j];
        k := 2 * (Length(A) - Rank(Concatenation(A, B)));
        if targetK > 0 and k <> targetK then
            continue;
        fi;
        encodingRate := k / n / 2;
        if encodingRate < encodingRateLowerBound then
            continue;
        fi;
        Hs := ConstructParityCheckMatrix(A, B);
        d := DistRandCSS(Hs[1], Hs[2], numInformationSets, distanceLowerBound, 2: field:=F);
        if d > 0 then
            Add(codes, rec(parameter := [n, k, d], As:=combs[i], Bs:=combs[j], encodingRate := encodingRate));
        fi;
    od;
    return codes;
end;


# Given the parameter l and m, search the bivariate bicyclic code
SearchBBCodes := function(l, m, resultsFilepath, encodingRateLowerBound, distanceLowerBound, numInformationSets, chunkSize, targetK)
    local x, y, group, gx, gy, sizeM, ps, pool, n, codes, allMats, combs, tasks, allPotentialGenerators;
    codes := [];
    x := CyclicMat(l, m, false);
    y := CyclicMat(l, m, true);
    n := 2 * l * m;
    # Cyclic group associated with the monomials
    group := ConstructGroup(l, m);
    sizeM := group[1];
    gx := group[2];
    gy := group[3];

    ps := GetSearchPool(l, m, gx, gy, x, y);
    combs := ps[1];
    allMats := ps[2];
    allPotentialGenerators := ps[3];
    pool := ps[4];
    Print("Search pool size: ", Length(pool), "\n");

    # CalculateFor(pool, combs, allMats, allPotentialGenerators,  sizeM, n, Rat(encodingRateLowerBound), numInformationSets, distanceLowerBound, targetK);
    tasks := List(SplitListIntoChunks(pool, chunkSize), chunk -> RunTask(CalculateFor, chunk, combs, allMats, allPotentialGenerators,  sizeM, n, Rat(encodingRateLowerBound), numInformationSets, distanceLowerBound, targetK));
    codes := Concatenation(List(tasks, task -> TaskResult(task)));
    PrintCSV(resultsFilepath, codes);
    return codes;
end;

# l := 12;
# m := 6;
# minDist := 10;
# rejectEncodingRate := 1 / 25;
# SearchBBCodes(l, m, rejectEncodingRate, minDist, 10, 1000, "results/codes.csv");

# l := 6;
# m := 6;
# distanceLowerBound := 5;
# encodingRateLowerBound := 1 / 13;
# targetK := 0;
# SearchBBCodes(l, m, "results/codes.csv", encodingRateLowerBound, distanceLowerBound, 10, 5000, targetK);
