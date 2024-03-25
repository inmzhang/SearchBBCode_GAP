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
    return [gx, gy];
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


FilterToricLayout := function(ijList, allPotentialGenerators, lm)
    local filteredList, ij;
    filteredList := [];
    for ij in ijList do
        if IsToricLayout(allPotentialGenerators[ij[1]], allPotentialGenerators[ij[2]], lm) then
            Add(filteredList, ij);
        fi;
    od;
    return filteredList;
end;


FilterEncodingRate := function(ijList, allMats, n, encodingRateLowerBound, targetK)
    local codes, ij, i, j, A, B, encodingRate, k;
    codes := [];
    for ij in ijList do
        i := ij[1];
        j := ij[2];
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
        Add(codes, rec(i:=i, j:=j, n:=n, k:=k, encodingRate := encodingRate));
    od;
    return codes;
end;


ComputeDistance := function(code, allMats, numInformationSets, distanceLowerBound)
    local Hs;
    Hs := ConstructParityCheckMatrix(allMats[code.i], allMats[code.j]);
    return DistRandCSS(Hs[1], Hs[2], numInformationSets, distanceLowerBound, 2: field:=F);
end;


# Given the parameter l and m, search the bivariate bicyclic code
SearchBBCodes := function(l, m, resultsFilepath, encodingRateLowerBound, distanceLowerBound, numInformationSets, chunkSize, targetK)
    local x, y, group, gx, gy, sizeM, ps, pool, n, codes, allMats, combs, tasks, allPotentialGenerators, satisfiedCodes, distances, d, code, i;
    codes := [];
    x := CyclicMat(l, m, false);
    y := CyclicMat(l, m, true);
    n := 2 * l * m;
    # Cyclic group associated with the monomials
    group := ConstructGroup(l, m);
    gx := group[1];
    gy := group[2];

    ps := GetSearchPool(l, m, gx, gy, x, y);
    combs := ps[1];
    allMats := ps[2];
    allPotentialGenerators := ps[3];
    pool := ps[4];
    Print("Search pool size: ", Length(pool), "\n");

    # Filter toric layout codes
    tasks := List(SplitListIntoChunks(pool, chunkSize), chunk -> RunTask(FilterToricLayout, chunk, allPotentialGenerators, l * m));
    codes := Concatenation(List(tasks, task -> TaskResult(task)));
    Print("Number of codes with toric layout: ", Length(codes), "\n");

    # Filter encoding rate
    tasks := List(SplitListIntoChunks(codes, Int(chunkSize / 2)), code -> RunTask(FilterEncodingRate, code, allMats, n, Rat(encodingRateLowerBound), targetK));
    codes := Concatenation(List(tasks, task -> TaskResult(task)));
    Print("Number of code satisfying encoding rate conditions: ", Length(codes), "\n");

    # Compute the distance
    tasks := List(codes, code -> RunTask(ComputeDistance, code, allMats, numInformationSets, distanceLowerBound));
    distances := List(tasks, task -> TaskResult(task));

    satisfiedCodes := [];
    for i in [1..Length(codes)] do
        code := codes[i];
        d := distances[i];
        if d > 0 then
            Add(satisfiedCodes, rec(parameter:=[code.n, code.k, d], encodingRate:=code.encodingRate, As:=combs[code.i], Bs:=combs[code.j]));
        fi;
    od;
    
    Print("Number of code satisfying all conditions: ", Length(satisfiedCodes), "\n");

    if Length(satisfiedCodes) > 0 then
        PrintCSV(resultsFilepath, satisfiedCodes);
    fi;
    return codes;
end;

# l := 12;
# m := 6;
# distanceLowerBound := 10;
# encodingRateLowerBound := 1 / 25;
# targetK := 0;
# SearchBBCodes(l, m, "out/codes.csv", encodingRateLowerBound, distanceLowerBound, 10, 5000, targetK);

# l := 6;
# m := 6;
# distanceLowerBound := 5;
# encodingRateLowerBound := 1 / 13;
# targetK := 0;
# SearchBBCodes(l, m, "out/codes.csv", encodingRateLowerBound, distanceLowerBound, 10, 5000, targetK);
