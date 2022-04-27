clear all
close all
clc
%%
ax=0.0;
bx=1.0;
Nx=20;
Nt=20;
number_mesh=4;
number_mesh_point=zeros(number_mesh,1);
norm_max=zeros(number_mesh,1);
norm_l2=zeros(number_mesh,1);
norm_maxh1=zeros(number_mesh,1);
norm_h1=zeros(number_mesh,1);
%% Solve discrite solution and refine mesh
for inumber_mesh=1:number_mesh
    number_mesh_point(inumber_mesh)=Nx;
    h=(bx-ax)/Nx;
    k=h^2; %compute k with the time step restriction  
    T=k*Nt; %compute T
%% Creat the mesh point
    x=zeros(Nx+1,1);
    for i=1:Nx+1
        x(i)=(i-1)*h;
    end
    
    t=zeros(Nt+1,1);
    for i=1:Nt+1
        t(i)=(i-1)*k;
    end
%% Creat matrix A
    r=1;
    A=sparse(Nx-1,Nx-1);
    for i=1:Nx-1
        if i==1
            A(i,i)=-2*r;
            A(i,i+1)=r;
        elseif i==Nx-1
            A(i,i)=-2*r;
            A(i,i-1)=r;
        else
            A(i,i-1)=r;
            A(i,i)=-2*r;
            A(i,i+1)=r;
        end
    end
%% Discrete 
    u0=zeros(Nx-1,1);
    for i=2:Nx
        u0(i-1)=u_exact(x(i),0);
    end
    F=zeros(Nx-1,1);
    u=zeros(Nx-1,1);
    for time=1:Nt+1
        %------------------------------------------------------------------
        % Creat matrix F
        %------------------------------------------------------------------
        for i=2:Nx
            F(i-1)=k*f(x(i),t(time));
        end
        %------------------------------------------------------------------
        % Discrete the soluiton
        %------------------------------------------------------------------
        u=inv(eye(Nx-1)-A)*(u0+F);
        u0=u;
    end
    %----------------------------------------------------------------------
    % Exact the solution
    %----------------------------------------------------------------------
    uex=zeros(Nx+1,1);
    for i=1:Nx+1
        uex(i)=u_exact(x(i),T);
    end
    %----------------------------------------------------------------------
    % Create discrete solution with boundary 
    %----------------------------------------------------------------------
    udis=zeros(Nx+1,1);
    for j=1:Nx+1
        if j==1
            udis(j)=0;
        elseif j==Nx+1
            udis(j)=0;
        else
            udis(j)=u(j-1);
        end
    end
%% Calculate the error on L^infinity

    norm_max(inumber_mesh)=0.0;
    for i=1:Nx+1
        if (abs(udis(i)-uex(i)) > norm_max(inumber_mesh))
            norm_max(inumber_mesh)=abs(udis(i)-uex(i));
        end
    end
    
    norm_max(inumber_mesh)
%%  Calculate the error on L^2 

    norm_l2(inumber_mesh)=0;
    for i=1:Nx+1
        norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+(udis(i)-uex(i))^2*h;
    end
    
    norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
    norm_l2(inumber_mesh)
%% Calculate the error on maxH1    

    norm_maxh1(inumber_mesh)=0;
    for i=1:Nx
        if (abs(((udis(i+1)-uex(i+1))-(udis(i)-uex(i)))/h) > norm_maxh1(inumber_mesh))
            norm_maxh1(inumber_mesh)=abs(((udis(i+1)-uex(i+1))-(udis(i)-uex(i)))/h);
        end
    end
    norm_maxh1(inumber_mesh)

%% Calculate the error on H1

    norm_h1(inumber_mesh)=0;
    for i=1:Nx
        norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+(((udis(i+1)-uex(i+1))-(udis(i)-uex(i)))/h)^2*h;
    end
    norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
%% Figure exact and dicrete solutions 
    figure
    plot(x, uex,'cyan h', x, udis,'magenta')
    legend('Exact Solution','Discrete Solution')
    xlabel('x');ylabel('value');
    title(['Comparison between exact and discrete solutions with Nx = ',num2str(Nx)]);
%% Refine mesh (increse mesh point)    
    Nx=Nx*2;
end

%% Figure for errors respect to number of mesh point
figure
plot(log(number_mesh_point), -log(norm_max),'blue', log(number_mesh_point), -log(norm_l2), 'red',...
    log(number_mesh_point), -log(norm_maxh1), 'cyan', log(number_mesh_point), -log(norm_h1), 'magenta',...
    log(number_mesh_point), 2*log(number_mesh_point),'black');
xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors');
legend('norm_max','norm_l2','norm_maxh1','norm_h1','2x','Location','NorthEastOutside');  


