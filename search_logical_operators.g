LoadPackage("QDistRnd");

F := GF(2);

codes := [
    [72, 12, 6, 6, 6],
    # [90, 8, 10, 15, 3],
    # [108, 8, 10, 9, 6],
    # [144, 12, 12, 12, 6],
    # [288, 12, 18, 12, 12],
];

loadParityCheckMats := function(hxFilePath, hzFilePath)
    local hxInfo, hzInfo, Hx, Hz;
    hxInfo := ReadMTXE(hxFilePath);
    hzInfo := ReadMTXE(hzFilePath);
    Hx := hxInfo[3];
    Hz := hzInfo[3];
    return [Hx, Hz];
end;


circularShiftVec := function(vec, s)
    local len, shift;
    len := Length(vec);
    shift := s mod len;
    return Concatenation(vec{[(len - shift + 1)..len]}, vec{[1..(len - shift)]});
end;


blockCircularShiftVec := function(vec, blockLength, innerShift, outerShift)
    local i, numBlock, innerShifted;
    numBlock := Int(Length(vec) / blockLength);
    innerShifted := Concatenation(List([1..numBlock], i -> circularShiftVec(vec{[(i-1)*blockLength+1..i*blockLength]}, innerShift)));
    return circularShiftVec(innerShifted, blockLength * outerShift);
end;


circularShiftedMat := function(vec, l, m)
    return Matrix(List(Cartesian([1..l], [1..m]), i -> blockCircularShiftVec(vec, m, i[2]-1, i[1]-1)));
end;


horizontalConcatenation := function(mat1, mat2)
    return Matrix(List([1..Length(mat1)], i -> Concatenation(mat1[i], mat2[i])));
end;


sparsePolynomialRepr := function(gf2Vec)
    local sparse, i;
    sparse := [];
    for i in [1..Length(gf2Vec)] do
        if gf2Vec[i] = One(F) then
            Add(sparse, i);
        fi;
    od;
    return sparse;
end;


sparseToVec := function(sparse, length)
    local vec, i;
    vec := ZeroVector(F, length);
    for i in sparse do
        vec[i] := One(F);
    od;
    return vec;
end;


getInversePolyIndex := function(i, l, m)
    local px, py, invpx, invpy;
    px := Int((i-1) / m);
    py := Int((i-1) mod m);
    invpx := ((-px) mod l) + 1;
    invpy := ((-py) mod m) + 1;
    return (invpx-1) * m + invpy;
end;


getInversePolyVec:= function(vec, l, m)
    local sparse;
    sparse := sparsePolynomialRepr(vec);
    return sparseToVec(List(sparse, i -> getInversePolyIndex(i, l, m)), Length(vec));
end;


polyVecProduct := function(a, b, m)
    local i, ret;
    ret := ZeroVector(F, Length(a));
    for i in [1..Length(a)] do
        if a[i] = Zero(F) then continue; fi; 
        ret := ret + blockCircularShiftVec(b, m, Int((i-1) mod m), Int((i-1) / m));
    od;
    return ret;
end;


checkNumLogicalQubits := function(H, f, g, h, k, l, m, isZ, zeroPadMat, rankH)
    local sf, sg, sh, invh, invf, invg, L1, L2;
    if isZ then
        invf := getInversePolyVec(f, l, m);
        invg := getInversePolyVec(g, l, m);
        invh := getInversePolyVec(h, l, m);
        L1 := horizontalConcatenation(circularShiftedMat(invh, l, m), circularShiftedMat(invg, l, m));
        L2 := horizontalConcatenation(zeroPadMat, circularShiftedMat(invf, l, m));
    else
        L1 := horizontalConcatenation(circularShiftedMat(f, l, m), zeroPadMat);
        L2 := horizontalConcatenation(circularShiftedMat(g, l, m), circularShiftedMat(h, l, m));
    fi;
    return Rank(Concatenation(L1, L2, H)) = k + rankH;
end;


searchLogicalOperators := function(Hx, Hz, n, k, d, l, m, searchNumLimit)
    local n2, rankHx, rankHz, foundList, A, B, fSpaceBasis, ghSpaceBasis, f, gh, g, h, PadNullMat, X1s, X2s, r;
    rankHx := Rank(Hx);
    rankHz := Rank(Hz);
    foundList := [];
    n2 := Int(n / 2);
    A := List(Hx, row -> row{[1..n2]});
    B := List(Hx, row -> row{[n2+1..n]});
    # X_a = X(af, 0), Z_a = Z(ah.T, ag.T), X_b = X(bg, bh), Z_b = Z(0, bf.T)
    # f, g, h satisfy Bf = 0, gB + hA = 0
    fSpaceBasis := Basis(VectorSpace(F, NullspaceMat(B)));
    Print("fSpaceBasis: ", Length(fSpaceBasis), "\n");
    ghSpaceBasis := Basis(VectorSpace(F, NullspaceMat(Concatenation(B, A))));
    Print("ghSpaceBasis: ", Length(ghSpaceBasis), "\n");
    PadNullMat := ZeroMatrix(F, n2, n2);
    for f in IteratorByBasis(fSpaceBasis) do
        if WeightOfVector(f) < d then continue; fi;
        for gh in IteratorByBasis(ghSpaceBasis) do
            if WeightOfVector(gh) < d then continue; fi;
            g := gh{[1..n2]};
            h := gh{[n2+1..n]};
            if WeightOfVector(g) = 0 or WeightOfVector(h) = 0 then continue; fi;
            # determine whether this set of operators span all k logical qubits
            if not checkNumLogicalQubits(Hx, f, g, h, k, l, m, false, PadNullMat, rankHx) then continue; fi;
            # This check should not be required due to symmetry
            # Anyway, it does not hurt to put it here
            if not checkNumLogicalQubits(Hz, f, g, h, k, l, m, true, PadNullMat, rankHz) then 
                Error("Hx satisfy but not Hz!");
            fi;
            Add(foundList, rec(n:= n, f:=sparsePolynomialRepr(f), g:=sparsePolynomialRepr(g), h:=sparsePolynomialRepr(h)));
            Print(foundList, "\n");
            if Length(foundList) = searchNumLimit then
                return foundList;
            fi;
        od;
    od;
    return foundList;
end;


searchMonomialsBasedOnOperator := function(f, h, l, m, k)
    local fh, length, ns, ms, n, mm, nfh, mt, fhmt;
    fh := polyVecProduct(f, h, m);
    length := Length(fh);
    ns := [];
    ms := [];
    for n in [1..length] do
        nfh := polyVecProduct(fh, sparseToVec([n], length), m);
        if ForAny(ms, i -> nfh[i] = One(F)) then continue; fi;
        Add(ns, n);
        for mm in [1..length] do
            if mm in ms then continue; fi;
            mt := getInversePolyIndex(mm, l, m);
            fhmt := polyVecProduct(fh, sparseToVec([mt], length), m);
            if fhmt[getInversePolyIndex(ns[Length(ns)], l, m)] = Zero(F) then continue; fi;
            if ForAny(ns{[1..Length(ns)-1]}, i -> fhmt[getInversePolyIndex(i, l, m)] = One(F)) then continue; fi;
            Add(ms, mm);
            break;
        od;
        if Length(ns) = Int(k / 2) then
            if Length(ms) <> Int(k / 2) then
                Error("Unreachable!");
            fi;
            return [ns, ms];
        fi;
    od;
end;


checkLogicalQubitsMonomials := function(ns, ms, f, h, l, m)
    local length, fh, i, j, nn, mm, nfh;
    fh := polyVecProduct(f, h, m);
    length := Length(fh);

    for i in [1..Length(ns)] do
        nn := ns[i];
        mm := ms[i];
        nfh := polyVecProduct(fh, sparseToVec([nn], length), m);
        if nfh[mm] = Zero(F) then return false; fi;
        for j in [1..Length(ms)] do
            if j = i then continue; fi;
            if nfh[ms[j]] = One(F) then return false; fi;
        od;
    od;
    return true;
end;


records := [];
searchNumLimit := 1;

for code in codes do
    n := code[1];
    k := code[2];
    d := code[3];
    l := code[4];
    m := code[5];
    Print("n = ", n, ", k = ", k, ", d = ", d, "\n");
    hxFilePath := Concatenation("bbcode_mtxe/", String(n), "_", String(k), "_", String(d), "_Hx", ".mtxe");
    hzFilePath := Concatenation("bbcode_mtxe/", String(n), "_", String(k), "_", String(d), "_Hz", ".mtxe");
    pms := loadParityCheckMats(hxFilePath, hzFilePath);

    startTime := Runtime();
    found := searchLogicalOperators(pms[1], pms[2], n, k, d, l, m, searchNumLimit);
    endTime := Runtime();
    elapsedTime := endTime - startTime;
    Print("Time taken: ~", StringTime(elapsedTime), "\n");

    csvFilePath := Concatenation("results/logical_operators_", String(n), "_", String(k), "_", String(d), ".csv");
    PrintCSV(csvFilePath, found);
    Append(records, found);
od;

PrintCSV("results/logical_operators_by_gap.csv", records);


# Test logical qubits check
# code := [144, 12, 12, 12, 6];
# n := code[1];
# k := code[2];
# d := code[3];
# l := code[4];
# m := code[5];

# n2 := Int(n / 2);
# PadNullMat := ZeroMatrix(F, n2, n2);

# f := sparseToVec([1, 7, 10, 13, 19, 34, 37, 43, 46, 49, 55, 70], n2);
# g := sparseToVec([3, 5, 7, 9, 14, 16], n2);
# h := sparseToVec([1, 2, 3, 4, 8, 10], n2);

# hxFilePath := Concatenation("bbcode_mtxe/", String(n), "_", String(k), "_", String(d), "_Hx", ".mtxe");
# hzFilePath := Concatenation("bbcode_mtxe/", String(n), "_", String(k), "_", String(d), "_Hz", ".mtxe");
# pms := loadParityCheckMats(hxFilePath, hzFilePath);

# rankHx := Rank(pms[1]);
# rankHz := Rank(pms[2]);
# Print(checkNumLogicalQubits(pms[1], f, g, h, k, l, m, false, PadNullMat, rankHx), "\n");
# Print(checkNumLogicalQubits(pms[2], f, g, h, k, l, m, true, PadNullMat, rankHz), "\n");


# # Search logical qubits monomials
# Print(searchMonomialsBasedOnOperator(f, h, l, m, k), "\n");

# # # Test logical qubits monomials
# ns := [1, 2, 14, 18, 21, 25];
# ms := [2, 6, 8, 1, 25, 33];
# Print(checkLogicalQubitsMonomials(ns, ms, f, h, l, m), "\n");
