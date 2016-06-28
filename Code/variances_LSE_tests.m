%signal_averaging; theoretic comparison: LS eastimator with and without covariance matrix.
x=linspace(0,10,50);
y=3*x+8*sin(x)+7*cos(x)+1*cos(5*x)+2*sin(3*x);

m=y+randn(size(y)).*2.5


R=[x;cos(x);sin(x);cos(5*x);sin(3*x)]';


beta=pinv(R)*m.'
ls1=beta.'*R.'

figure(1)
hold on
plot(x,m,'+')
plot(x,y)
plot(x,ls1,'k-')
hold off
%%

m2=y+randn(size(y)).*(1.25+x./4)


beta=pinv(R)*m2.'
ls2=beta.'*R.'



figure(2)
hold on
plot(x,m2,'+')
plot(x,y)
plot(x,ls2,'k-')

hold off