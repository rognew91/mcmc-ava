function [pa1, pa2, pa3, pa4, pa5, pa6, pa7, pa8, pa9, pa10] = plot_reference(m_best, mode_model, medians, iqr_model, scale)
%function plot_reference(m_best, medians, iqr_model, scale)

%PLOT_REFERENCE Summary of this function goes here
%   Detailed explanation goes here

load('reference.mat')

%% make a unit sphere to base everything else on
[x,y,z] = sphere(30);

%% ice

p1 = scale*x + reference.alpha2(1);
q1 = scale*y + reference.beta2(1);
r1 = scale*z + reference.rho2(1);

fvc1 = surf2patch(p1,q1,r1);
pa1 = patch(fvc1, 'FaceColor', [0 0 1], 'EdgeColor','none');
reducepatch(pa1,0.5)
%s1 = surf(p1,q1,r1);

%% basement

p2 = scale*x + reference.alpha2(2);
q2 = scale*y + reference.beta2(2);
r2 = scale*z + reference.rho2(2);

%s2 = surf(p2,q2,r2);
fvc2 = surf2patch(p2,q2,r2);
pa2= patch(fvc2, 'FaceColor', [1 0 0], 'EdgeColor','none');
reducepatch(pa2,0.5)


%% lithified seds

p3 = scale*x + reference.alpha2(3);
q3 = scale*y + reference.beta2(3);
r3 = scale*z + reference.rho2(3);

%s3 = surf(p3,q3,r3);
fvc3 = surf2patch(p3,q3,r3);
pa3=patch(fvc3, 'FaceColor', 'm', 'EdgeColor','none');
reducepatch(pa3,0.5)


%% dilatant seds

p4 = scale*x + reference.alpha2(4);
q4 = scale*y + reference.beta2(4);
r4 = scale*z + reference.rho2(4);

%s4 = surf(p4,q4,r4);
fvc4 = surf2patch(p4,q4,r4);
pa4=patch(fvc4, 'FaceColor', 'g', 'EdgeColor','none');
reducepatch(pa4,0.5)


%% water

p5 = scale*x + reference.alpha2(5);
q5 = scale*y + reference.beta2(5);
r5 = scale*z + reference.rho2(5);

%s5 = surf(p5,q5,r5);
fvc5 = surf2patch(p5,q5,r5);
pa5=patch(fvc5, 'FaceColor', 'b', 'EdgeColor','none');
reducepatch(pa5,0.5)


%% permafrost

p6 = scale*x + reference.alpha2(6);
q6 = scale*y + reference.beta2(6);
r6 = scale*z + reference.rho2(6);

%s6 = surf(p6,q6,r6);
fvc6 = surf2patch(p6,q6,r6);
pa6=patch(fvc6, 'FaceColor', 'y', 'EdgeColor','none');
reducepatch(pa6,0.5)


%% stiff till

p7 = scale*x + 1800;
q7 = scale*y + 1000;
r7 = scale*z + 1900;

%s7 = surf(p7,q7,r7);
fvc7 = surf2patch(p7,q7,r7);
pa7 = patch(fvc7, 'FaceColor', [0.8500 0.3250 0.0980], 'EdgeColor','none');
reducepatch(pa7,0.5)


%% best model

p8 = scale*x + mode_model(4);
q8 = scale*y + mode_model(6);
r8 = scale*z + mode_model(2);

%s8 = surf(p8, q8, r8);
fvc8 = surf2patch(p8,q8,r8);
pa8=patch(fvc8, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor','none');
reducepatch(pa8,0.5)


%% median
if ~isnan(sum(medians)) || ~isnan(sum(iqr_model))

    p9 = 0.5*iqr_model(4)*x + medians(4);
    q9 = 0.5*iqr_model(6)*y + medians(6);
    r9 = 0.5*iqr_model(2)*z + medians(2);

    %s9 = surf(p9, q9, r9);
    fvc9 = surf2patch(p9,q9,r9);
    pa9=patch(fvc9, 'FaceColor', 'w', 'FaceAlpha', 0.8, 'EdgeColor', 'none');
    reducepatch(pa9,0.5)


else 
    pa9 = NaN;

end

%shading interp
camlight
%lighting phong


    p10 = scale*x + m_best(4);
    q10 = scale*y + m_best(6);
    r10 = scale*z + m_best(2);

    %s9 = surf(p9, q9, r9);
    fvc10 = surf2patch(p10,q10,r10);
    pa10=patch(fvc10, 'FaceColor', 'k', 'EdgeColor', 'none');
    reducepatch(pa10,0.5)
end
