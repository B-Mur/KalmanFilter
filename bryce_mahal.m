function dist  = bryce_mahal(a, A, b, B)

dist = sqrt((a-b)' * inv(A+B) * (a-b));

