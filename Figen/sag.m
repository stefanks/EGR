%SAG with quadratic f

m = 100; %component functions
n = 20; %dimension

%create random component functions
H = zeros(n,n,m);
q = zeros(n,m);
barH = zeros(n);
for i=1:m
    h = rand(n);
    H(:,:,i) = i*h'*h;
    q(:,i) = rand(n,1);
    barH = barH + H(:,:,i);
end
barq = sum(q,2)/m;
barH = barH/m;

L = norm(barH);
mu = min(eig(barH));

x_true = barH\barq;

%var_g_x_true = 0;
%for i=1:m
%    var_g_x_true = var_g_x_true+(H(:,:,i)*x_true-q(:,i))'*(H(:,:,i)*x_true-q(:,i));
%end
%var_g_x_true = var_g_x_true/m;

x = m*rand(n,1); %initial point

y = zeros(n,m); %storage of gradient evaluations
for i=1:m
    y(:,i) = H(:,:,i)*x - q(:,i);
end

maxiter = 1000;
e = zeros(maxiter,1);
dx = zeros(maxiter,1);
nv = zeros(maxiter,1);
Mv = zeros(maxiter,1);
f = zeros(maxiter,1);

eta = 0.1;
%alpha = eta / L;
alpha = 0.25* eta / L;

x = x - alpha*sum(y,2)/m;

for k=1:maxiter-1
    
    %update y
    %selected = 0;
    %while(selected<eta*m)
    for i=1:m
        if(rand<eta)
            %i = randint(1,1,[1,m]);
            y(:,i) = H(:,:,i)*x - q(:,i);
            %selected = selected+1;
        end
    end
    
    x = x - alpha*sum(y,2)/m;
            
    %update errors
    e(k) = norm(sum(y,2)/m-(barH*x-barq));
    dx(k) = norm(x - x_true);
        
    M = [(1-eta)*(eye(n)+alpha*barH) , (1-eta)*alpha*barH*barH/L;  
                 -alpha*eye(n)*L         ,    eye(n)-alpha*barH];

    v = [(sum(y,2)/m-(barH*x-barq))/L;(x - x_true)];
    nv(k) = norm(v);
    Mv(k) = norm(M*v);
    
    if(Mv(k)>nv(k))
        disp('Increase in v!');
    end
    
    %f(k)=0;
    %for i=1:m
    %    f(k)=f(k)+1/2*x'*H(:,:,i)*x - x'*q(:,i);
    %end
    %f(k) = f(k)/m;   
    
    fprintf('iter:%d, norm(Mv):%f, norm(v):%f, dx:%f, e:%f\n',k,Mv(k),nv(k),dx(k),e(k));
        
end
%norm(sqrt(M'*M))
plot(log(nv),'r'); hold all
plot(log(Mv),'b');
%plot(log(f),'k');
plot(log(dx),'g');