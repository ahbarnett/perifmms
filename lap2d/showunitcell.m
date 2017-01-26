function showunitcell(U)
e1=U.e1; e2 = U.e2;
plot(-(e1+e2)/2 + [0 e1 e1+e2 e2 0],'k-');
