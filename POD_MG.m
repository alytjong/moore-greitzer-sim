load Yfig7.mat;

%Yfig6.mat is stable
%Yfig7.mat is surge
%Yfig8.mat is stall
%Yfig9.mat is surge + stall
%%
L = 2*pi;
n = 512;  %domain, frequency resolution

tht2 = linspace(-L/2, L/2, n+1);
tht = tht2(1:n);%Check mode shape

t = 0:1:1000;

%%
gtsol = zeros(numel(t), numel(tht));
for j = 1:length(t)
    gtsol(j,:) = real(ifft(Y(j,1:end-2)));
end
%%
% figure;
% line1 = plot(tht,gtsol(1,:));
% for iii = 1:size(gtsol,1)
%     set(line1, 'ydata', gtsol(iii,:));
%     axis([-pi pi -10 10]);
%     title(['g(\theta, \xi=' num2str(t(iii)) ')']);
%     drawnow;
%     pause(0.01);
% end

%%  Perform SVD
G = gtsol';
%SVD
[U,S,V] = svd(G);

%Check for dominant modes
figure('color', 'w');
pcent = diag(S)*100./sum(diag(S));
bar(pcent(1:10));
title('Percentage of top 10 energy modes of the stable PDE solution', ...
    'interpreter', 'latex', 'fontsize', 18);

%%
%Found r dominant modes
r = 16;
Phi = U(:,1:r);

% %Verify orthonormality
% Phi'*Phi

%Calculate amplitudes
%Convert to double loop
a = zeros(r,size(G,2));
for iii = 1:size(G,2)
    for jjj = 1:r
        %a1 = dot(G(:,iii),Phi(:,1))/dot(Phi(:,1),Phi(:,1));
        %a2 = dot(G(:,iii),Phi(:,2))/dot(Phi(:,2),Phi(:,2));
        a(jjj,iii) = dot(G(:,iii),Phi(:,jjj))/dot(Phi(:,jjj),Phi(:,jjj));
    end
    %a(:,iii) = [a1; a2];
end

%%  Results

%Check amplitudes over time

% %Stable figure
% figure('color', 'w');
% for iii = 1:r
%     subplot(r,1,iii); plot(t, a(iii,:));
%     if iii == 1
%         title('Reduced order amplitudes of the stable case, $$r=2$$', ...
%     'interpreter', 'latex', 'fontsize', 18);
%     end
%     ylabel(['$$a_' num2str(iii) '(\xi)$$'], 'interpreter', 'latex', ...
%         'fontsize', 15);
% end
% xlabel('Time $$\xi$$', 'interpreter', 'latex', 'fontsize', 15);


%Surge, stall, combination figures
figure('color', 'w', 'position', [360 70 490 636]);
for iii = 1:r
    subplot(r,1,iii); plot(t(1:300), a(iii,1:300));
    if iii == 1
        title('Reduced order amplitudes of the combination case, $$r=16$$', ...
    'interpreter', 'latex', 'fontsize', 18);
    end
    ylabel(['$$a_{' num2str(iii) '}(\xi)$$'], 'interpreter', 'latex', ...
        'fontsize', 10);
    if iii ~= r
        set(gca,'xticklabel',[])
    end
    set(gca, 'yticklabel', []);
end
xlabel('Time $$\xi$$', 'interpreter', 'latex', 'fontsize', 15);

%%
%Check mode shapes over time

% %Stable case
% figure('color', 'w');
% for iii = 1:r
%     subplot(r,1,iii); plot(tht, Phi(:,iii));
%     if iii == 1
%         title('Reduced order mode shapes of the stable case, $$r=2$$', ...
%     'interpreter', 'latex', 'fontsize', 18);
%     end
%     ylabel(['$$U_' num2str(iii) '(\theta)$$'], 'interpreter', 'latex', ...
%         'fontsize', 10);
%     xticks([-pi 0 pi]);
%     xticklabels({'-\pi','0','\pi'});
% end
% xlabel('$$\theta$$', 'interpreter', 'latex', 'fontsize', 15);

%Surge figure
figure('color', 'w', 'position', [360 70 490 636]);
for iii = 1:r
    subplot(r,1,iii); plot(tht, Phi(:,iii));
    if iii == 1
        title('Reduced order mode shapes of the combination case, $$r=16$$', ...
    'interpreter', 'latex', 'fontsize', 18);
    end
    ylabel(['$$U_{' num2str(iii) '}(\theta)$$'], 'interpreter', 'latex', ...
         'fontsize', 10);
    if iii ~= r
        set(gca,'xticklabel',[])
    end
    set(gca, 'yticklabel', []);
    if iii == r
        xticks([-pi 0 pi]);
        xticklabels({'-\pi','0','\pi'});
    end
end
xlabel('$$\theta$$', 'interpreter', 'latex', 'fontsize', 15);



%%
%Check POD result with real data
G_approx = zeros(size(G));
for iii = 1:numel(t)
    a_tht = zeros(numel(tht),1);
    for jjj = 1:r
        a_tht = a_tht + a(jjj,iii)*Phi(:,jjj);
    end
    G_approx(:,iii) = a_tht;
end

f = figure('color', 'w'); hold on;
line1 = plot(tht,G(:,1),'b');
line2 = plot(tht, G_approx(:,1),'r.');
title('Dynamics of $$g(\xi,\theta)$$ in surge configuration', ...
    'interpreter', 'latex', 'fontsize', 18);
xlabel('$$\theta$$' ,'interpreter', 'latex', 'fontsize', 15);
ylabel('$$g(\xi,\theta)$$', 'interpreter', 'latex', 'fontsize', 15);
xticks([-pi 0 pi]);
xticklabels({'-\pi','0','\pi'});
for iii = 1:numel(t)
    set(line1, 'ydata', G(:,iii));
    set(line2, 'ydata', G_approx(:,iii));
    axis([-pi pi -3 3]);
    legend('Original data', 'Reduced data');
    drawnow;
    pause(0.01);
    
    filename = ['podanimfig7_' num2str(iii) '.png'];
    saveas(f, filename);
end






