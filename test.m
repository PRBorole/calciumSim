A = [1 2 3 4; 5 6 7 8; 9 10 11 12; 13 14 15 16];

n = 60000;
B = zeros(n*4);
for i = 1:4:n*4
    B(i:i+3,i:i+3) = A;
end