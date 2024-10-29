%This program was coded by Erfan Radfar
%precision is the minimum step of Omega
%max_iteration is the maximum number of points in the plot
function Output=nyquistfull(G,precision,max_iteration)
syms s;
step=precision*2; %initializing of the first step size of Omega
counter=1;
max_step=precision^(-2);%maximum allowable step
data=[];          %complex points of plot G(iw) are stored here
[Num,Den] = tfdata(G);
system = poly2sym(cell2mat(Num),s)/poly2sym(cell2mat(Den),s);%system is simbolic expression of G
while (step>precision && step<max_step)
    
    if counter==1
        w=precision;
        data(counter,1)=double(subs(system,1i*w)) ;
        data(counter,2)= w;
    else
        w=step/2+w;
        data(counter,1)=double(subs(system,1i*w)) ;
        data(counter,2)= w;
        step = step *(  1/abs( data(counter)- data(counter-1) ) + 1)^(1/20);
        %now we update our step
        % 1 is added to denominator to stop the step from becoming too small
        % in case of so small difference in subsequent data
        
        
    end
    
    counter=counter+1;
    
    if counter>max_iteration
        break
    end
    
end

my_plot=plot(data(:,1),'b-o','MarkerSize',3);
my_plot.DataTipTemplate.DataTipRows(1).Label = "Real";
my_plot.DataTipTemplate.DataTipRows(2).Label = "Imaginary"; 
my_plot.DataTipTemplate.DataTipRows(3)=dataTipTextRow("Frequency(rad/s)",data(:,2));
grid on
hold on
xlabel('Im (\omega)');
ylabel('Re (\sigma)');
title('Nyquist Diagram');
my_plot=plot(conj(data(:,1)),'r-o','MarkerSize',3);
my_plot.DataTipTemplate.DataTipRows(1).Label = "Real";
my_plot.DataTipTemplate.DataTipRows(2).Label = "Imaginary"; 
my_plot.DataTipTemplate.DataTipRows(3)=dataTipTextRow("Frequency(rad/s)",-data(:,2));

hold on
plot_poles(G,data,precision,system); %This function, plots circles and semicircles
%which are for poles of G on imaginary axis
%True radiuses are infinit


Output=cat( 1 , flip( [conj(data(:,1)),-data(:,2)] ) , data);
%Column 1 is G(jw) and Column 2 is corresponding w
end




function abs_poles=overlap_check(zeros_,poles)
z=zeros_; p=poles;
m=length(z); n=length(p);L=1;
while L<=m
    k=1;
    while k<=n
        if abs(z(L)-p(k))<0.001
            p(k)=[];
            z(L)=[];
            m=m-1;
            k=k-1;
            n=n-1;
            
        end
        k=k+1;
    end
    L=L+1;
end


N=length(p);
temp=ones(1,N);%Used to store iteration of each pole
n=1;
while n<=N
    m=n;
    while m<=N
        
        if (abs(p(n)-p(m))<0.001) && (n~=m)
            temp(n)=temp(n)+1; 
            temp(m)=[];
            p(m)=[];
            N=N-1;m=m-1;
        end
        m=m+1;
        
    end
    
    n=n+1;
end
abs_poles=cat( 2 , p , transpose(temp) );
%Now we have every poles and thier iterations
end



function plot_poles(G,data,precision,system)
poles=pole(G);
zeros=zero(G);
abs_poles=overlap_check(zeros,poles);
%returns all the absolute poles of system
%by removing same roots at numerator and denominator
%Now we find the imaginary poles
[M,M_]=size(abs_poles);
im_poles=[];
for m=1:M
    if abs(real(abs_poles(m,1)))<0.001
       im_poles=cat( 1 , im_poles , abs_poles(m,:) );
    end
end
[H,H_]=size(im_poles);
nearest_omega=[];
[N,N_]=size(data);
for h=1:H
    min_distance=100;
    for n=1:N
        distance=abs(data(n,2)-abs_poles(h,1));
        if distance < min_distance
            % checks and find the the omega
            % which is nearest to the pole and saves it in nearest omega
            %then it will be used
            %for start and end point
            %of circles and
            %semicircles
            min_distance=distance;
            nearest_omega(h)=n;
        end
    end
end

%Now we find the furthest point from the origin  and call its distance R
R=max(abs(data(nearest_omega)));

theta=transpose((-pi/2) :0.01: (pi/2) );
for k=1:H
    w(:,k)=(precision*((theta-pi)/10+1)) .* exp(1i*theta) + im_poles(k,1);
    Circle_(:,k)=[double(subs(system,w(:,k))) ; data(nearest_omega(k),1)];
    Circle(:,k)=[ conj(data(nearest_omega(k),1)) ; Circle_(:,k) ];
    my_plot=plot(Circle(:,k),'g--o','MarkerSize',3);
    my_plot.DataTipTemplate.DataTipRows(1).Label = "Real";
    my_plot.DataTipTemplate.DataTipRows(2).Label = "Imaginary"; 
end

end
