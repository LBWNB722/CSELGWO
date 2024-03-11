function [gbestfitness,gbest,Convergence_curve]= CSELGWO(n,maxiter,lb,ub,dim,fobj)
Alpha_pos = zeros(1, dim);
Alpha_score = inf;
Beta_pos = zeros(1, dim);
Beta_score = inf;
Delta_pos = zeros(1, dim);
Delta_score = inf;
CR=0.9;
new_fitness=zeros(1,n);
archive_pos=zeros(n,dim);
levy_pos=zeros(n,dim);
X1=[];
X2=[];
X3=[];
ns=30;
nd=5;
jishu=0;
op=[];
pos = initialization(n,dim,ub,lb);
Convergence_curve = zeros(1, maxiter);
for i=1:n
    fitness(i)=fobj(pos(i,:));
end
[sortfitness,indexsort]=sort(fitness);
Alpha_pos=pos(indexsort(1),:);
Alpha_score=sortfitness(1);
Beta_pos=pos(indexsort(2),:);
Beta_score=sortfitness(2);
Delta_pos=pos(indexsort(3),:);
Delta_score=sortfitness(3);
gbest=Alpha_pos;
gbestfitness=Alpha_score;
archive_fitness=fitness;
archive_pos=pos;
t=1;
while t<maxiter
    a = 2 - t * 2 / maxiter;
    sizepop=size(pos(:,1));
    v=zeros(sizepop(1),dim);
    u=zeros(sizepop(1),dim);
    for i=1:sizepop(1)
        Flag4ub = pos(i, :) > ub;
        Flag4lb = pos(i, :) < lb;
        pos(i, :) = (pos(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
        for j=1:dim
            r1 = rand();
            r2 = rand();
            A1 = 2*a*r1 - a;
            C1 = 2*r2;
            D_alpha = abs(C1*Alpha_pos(j) - pos(i, j));
            X1(j) = (Alpha_pos(j)  - A1*D_alpha);
            r1 = rand();
            r2 = rand();
            A2 = 2*a*r1 - a;
            C2 = 2*r2;
            D_beta = abs(C2*Beta_pos(j) - pos(i, j));
            X2(j) = (Beta_pos(j) - A2*D_beta);
            r1 = rand();
            r2 = rand();
            A3 = 2*a*r1 - a;
            C3 = 2*r2;
            D_delta = abs(C3*Delta_pos(j) - pos(i, j));
            X3(j) =(Delta_pos(j)  - A3*D_delta);
        end
        pos(i, :) = ((X1 + X2 + X3) / 3);
        if fitness(i)<Alpha_score
            Alpha_score=fitness(i);
            Alpha_pos=pos(i,:);
        end
        if fitness(i) > Alpha_score && fitness(i) < Beta_score
            Beta_score = fitness(i);
            Beta_pos = pos(i, :);
        end
        if fitness(i) > Alpha_score && fitness(i) > Beta_score && fitness(i) < Delta_score
            Delta_score = fitness(i);
            Delta_pos = pos(i, :);
        end

        fitness(i)=fobj(pos(i,:));
    end
    for i=1:sizepop(1)
        if archive_fitness(i)<fitness(i)
            pos(i,:)=archive_pos(i,:);
            fitness(i)=archive_fitness(i);
        end
    end
    F(t)=1-abs(normrnd(n*atan(t)/maxiter,0.01,[1 1])-1/2);
    [m, ~] = size(pos);
    random_matrix = rand(m, n);
    [~, indices] = sort(random_matrix(:,sizepop(1)));
    for i=1:sizepop(1)
        rand_pos(i,:)=pos(indices(i),:);
        rand_fitness(i)=fitness(indices(i));
    end
    app(t) = ceil(ns - (ns - nd) * (t / maxiter));
    now_app=app(t);
    [pos_1,pos_2,pos_3,pos_4,jishu] = Subpop(sizepop(1),dim,pos,rand_fitness,now_app,jishu,t,maxiter);
    for i=1:sizepop(1)
        r1=randi([1,sizepop(1)],1,1);
        while(r1==i)
            r1=randi([1,sizepop(1)],1,1);
        end
        new_r1=randi([1,floor(sizepop(1)/4)],1,1);
        while(new_r1==i)
            new_r1=randi([1,floor(sizepop(1)/4)],1,1);
        end
        new_r2=randi([1,floor(sizepop(1)/4)],1,1);
        while(new_r2==new_r1)||(new_r2==i)
            new_r2=randi([1,floor(sizepop(1)/4)],1,1);
        end
        new_r3=randi([1,floor(sizepop(1)/4)],1,1);
        while(new_r3==i)||(new_r3==new_r2)||(new_r3==new_r1)
            new_r3=randi([1,floor(sizepop(1)/4)],1,1);
        end
        new_r4=randi([1,floor(sizepop(1)/4)],1,1);
        while(new_r4==i)||(new_r4==new_r2)||(new_r4==new_r1)||(new_r4==new_r3)
            new_r4=randi([1,floor(sizepop(1)/4)],1,1);
        end
        h=rand;
        R=0.5;
        if h>R
            v(i,:)=pos(r1,:)+F(t)*(pos_1(new_r1,:) -pos_2(new_r2,:));
        else
            v(i,:)=pos(r1,:)+F(t)*(pos_3(new_r3,:) -pos_4(new_r4,:));
        end
    end
    r=randi([1,dim],1,1);
    for j=1:dim
        cr=rand;
        if (cr<=CR)||(n==r)
            u(:,j)=v(:,j);
        else
            u(:,j)=pos(:,j);
        end
    end
    for i=1:sizepop(1)
        new_fitness(i)=fobj(u(i,:));
        if new_fitness(i)<fitness(i)
            pos(i,:)=u(i,:);
            fitness(i)=new_fitness(i);
        end
        if fitness(i)<Alpha_score
            Alpha_score=fitness(i);
            Alpha_pos=pos(i,:);
        end
    end
    for i=1:sizepop(1)
        Flag4ub = pos(i, :) > ub;
        Flag4lb = pos(i, :) < lb;
        pos(i, :) = (pos(i, :) .* (~(Flag4ub + Flag4lb))) + ub .* Flag4ub + lb .* Flag4lb;
    end
    is(t)=(maxiter-t+1)/maxiter;
    deta(t)=is(t)*(t/maxiter);
    if t>=2
        average_fitness=sum(Convergence_curve)/t;
        op(t)=(Convergence_curve(t-1)-(Alpha_score))/(abs(average_fitness-Alpha_score));
        if 0<=op(t) && op(t)<=deta(t)
            for i=1:sizepop(1)
                step=Levy(dim);
                levy_pos(i,:)=pos(i,:)+step.*(Alpha_pos-pos(i,:));
                levy_fitness(i)=fobj(levy_pos(i,:));
                if levy_fitness(i)<fitness(i)
                    pos(i,:)=levy_pos(i,:);
                    fitness(i)=levy_fitness(i);
                end
                if fitness(i)<Alpha_score
                    Alpha_score=fitness(i);
                    Alpha_pos=pos(i,:);
                end
            end
        end
    end
    archive_fitness=fitness;
    archive_pos=pos;
    gbest=Alpha_pos;
    gbestfitness=Alpha_score;
    Convergence_curve(t)=gbestfitness;
    t=t+1;
end
end





