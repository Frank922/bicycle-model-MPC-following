syms f(x) v eta lr lf
f(x) = v*cos(eta+atan((lr/(lf + lr))*tan(x)));
df_1 = diff(f,x);
f(x) = v*sin(eta+atan((lr/(lf + lr))*tan(x)));
df_2 = diff(f,x);
f(x) = (v/lr)*sin(atan((lr/(lf + lr))*tan(x)));
df_3 = diff(f,x);
beta = atan((lr/(lf+lr))*tan(x));
A = [0 0 -v*sin(eta+beta) cos(eta+beta); 0 0 v*cos(eta+beta) sin(eta+beta); 0 0 0 sin(beta)/lr; 0 0 0 0];
B = [df_1 0; df_2 0; df_3 0; 0 1];
