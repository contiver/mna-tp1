function comparison = compareEigenvalueMethods(m)
    [withWholeMatrix, _] = sort(eig(formIsingMatrix(m)), 'ascend');
    [with2by2Matrixes, _] = sort(getEigenvalues(m), 'ascend');
    comparison = [withWholeMatrix, with2by2Matrixes];
end

function eigenvalues = getEigenvalues(m)
    E = formE();
    eigenvalues = [];
    for k = 1:m
        F = formF(pi/4, m, k);
        K = E * F;
        [first, second] = get2by2MatrixEigenvalues(K);
        eigenvalues = [eigenvalues; first; second];
    end
end

function E = formE(alpha=pi/4)
    E = [cos(alpha), sin(alpha); -sin(alpha), cos(alpha)];
end

function F = formF(beta=pi/4, m, k)
    angle = formAngle(m);
    cosB = cos(beta);
    sinB = sin(beta);
    F = [cosB, -(angle ** k) * sinB; (angle ** (m - k)) * sinB, cosB];
end

function angle = formAngle(m)
    angle = e ** (i * 2 * pi / m);
end

function [first, second] = get2by2MatrixEigenvalues(A)
    traza = trace(A);
    discriminante = sqrt(traza ** 2 - 4 * det(A));
    first = (traza + discriminante) / 2;
    second = (traza - discriminante) / 2;
    [first, second];
end

function matrix = formIsingMatrix(m, alpha = pi / 4, beta = pi / 4)
    a = alpha;
    b = beta;
    E = [cos(a), sin(a); -sin(a), cos(a)];
    F = [cos(b), sin(b); -sin(b), cos(b)];

    K = createMatrix_K(m, E);
    L = createMatrix_L(m, F);
    matrix = K * L;
end

function matrix = createMatrix_K(m, E)
    matrix = insertMatrix(m, E);
end

function matrix = createMatrix_L(m, F)
    buffer = insertMatrix(m - 1, F);
    k = 2 * m;
    matrix = zeros(k);
    matrix(1,1) = matrix(k, k) = F(1,1);
    matrix(k, 1) = F(1, 2);
    matrix(1, k) = F(2, 1);
    matrix(2:(k - 1), 2:(k - 1)) = buffer(:,:);
end

function matrix = insertMatrix(m, E)
    matrix = zeros(2 * m);
    for i = 1:2:(2 * m)
        matrix(i:i+1, i:i+1) = E(:,:);
    end
end
