%1matrixid, 
%2 nodes,
%3.edges, non zero elements,
%4. gem fillins (lexp)
%5. gem processed rows (lexp)
%6. gem fillins (lexm)
%7. gem processed rows (lexm)
%8. gem fillins (orig)
%9. gem processed rows (orig)
%10. input space 
%11. space used by gem (lexp)
%12. space used by gem (lexm)
%13. space used by gem (orig)
%14, gem runtime(lexp)
%15, gem runtime(lexm)
%16, gem runtime(orig)
%17. time taken to initialize matrix
%close all;
data = [68, 2003, 11973, 0, 3, 0, 1, 0, 3, 32048, 16080, 16080, 16080, 0, 0, 0, 0.11;
    220, 100, 347, 986, 100, 938, 100, 916, 100, 1600, 41516, 39596, 38716, 0.01, 0 , 0, 0;
    240, 1000, 2375,  0, 311, 0, 305, 45220, 988, 16000, 8056, 8056, 1816856, 0, 0, 0.22, 0.01;
    344, 588, 12429, 59700, 588, 55710, 588, 53962, 588, 9408, 2399836, 2240236, 2170316, 0.76, 0.66, 0.62, 0.1;
    876, 306, 1162, 11728, 306, 11652, 306, 21878, 306, 4896, 475316, 472276, 881316, 0.05, 0.06, 0.16, 0.01;
    889, 9801, 48413, 2461474, 9801, 2461474, 9801, 1863176, 9801, 156816, 98655056, 98655056, 74723136, 45.24, 44.93, 22.89, 0.29;
    1205, 11445, 93781, 0, 0, 0, 0, 0, 0, 183120, 91616, 91616, 91616,0, 0, 0, 1.44;
    1239 3002 5000 0 1 0 0 0 0 48032 24096 24072 24072, 0, 0.01, 0, 0.06;
    1427 20000 30000 0 20000 0 20000 0 20000 320000 520076 520076 520076, 0.04, 0.03, 0.04, 0.19;
    1546, 4563, 17969, 33732 2313 7178 1658 1092580 4563  73008 1385840 323680 43794536 0.23 0.01 78.81 0.25;
    1553 6611 35472 70866 3151 13872 2549 2624446 6611 105776 2887584 607824 105110136 0.37 0.04 417.25 0.47;
    2401 1589 2742 0 0 0 0 0 0 25424 12768 12768 12768 0 0 0 0.02;
   ];
    
data(:,13)
inputmemory = data(:, 10); edges = data(:, 3); nodes = data(:, 2);
loglog(edges+nodes , inputmemory, 'o'); grid on;

gemmemory_lexp = data(:, 11); gemmemory_lexm = data(:, 12); gemmemory_orig = data(:, 13);
fillins_gem_lexp = data(:, 4);fillins_gem_lexm = data(:, 6);fillins_gem_orig = data(:, 8);
figure;
loglog(fillins_gem_lexp, gemmemory_lexp, 'ro', 'MarkerSize', 10);hold on;
loglog(fillins_gem_lexm, gemmemory_lexm, 'kx', 'MarkerSize', 10);hold on;
loglog(fillins_gem_orig, gemmemory_orig, 'b*', 'MarkerSize', 10);grid on;

inputtime = data(:, 17); edgesperrow = data(:,3)./data(:,2);
figure;
loglog(edges .* edgesperrow, inputtime, 'o'); grid on;

processed_lexp = data(:, 5); processed_lexm = data(:, 7); processed_orig = data(:, 9);
fillins_gem_lexp_perrow = fillins_gem_lexp./nodes;
fillins_gem_lexm_perrow = fillins_gem_lexm./nodes;
fillins_gem_orig_perrow = fillins_gem_orig./nodes;
gemtime_lexp = data(:, 14); gemtime_lexm = data(:, 15); gemtime_orig = data(:, 16);

finished_lexp = data(:, 5) == data(:, 2); finished_lexm = data(:, 7) == data(:, 2); finished_orig = data(:, 9) == data(:, 2);
solntie_lexp = finished_lexp .* nodes .* (fillins_gem_lexp_perrow + edgesperrow);
solntie_lexm = finished_lexm .* nodes .* (fillins_gem_lexm_perrow + edgesperrow);
solntie_orig = finished_orig .* nodes .* (fillins_gem_orig_perrow + edgesperrow);

%figure;
%loglog(processed_lexp .* (edgesperrow + fillins_gem_lexp_perrow), gemtime_lexp, 'o'); hold on;
%loglog(processed_lexm .* (edgesperrow + fillins_gem_lexm_perrow), gemtime_lexm, 'o'); hold on;
%loglog(processed_orig .* (edgesperrow + fillins_gem_orig_perrow), gemtime_orig, 'o'); hold on; grid on;

figure;
loglog(processed_lexp .* (edgesperrow + fillins_gem_lexp_perrow) + solntie_lexp, gemtime_lexp, 'ro', 'MarkerSize', 10); hold on;
loglog(processed_lexm .* (edgesperrow + fillins_gem_lexm_perrow) + solntie_lexm, gemtime_lexm, 'kx', 'MarkerSize', 10); hold on;
loglog(processed_orig .* (edgesperrow + fillins_gem_orig_perrow) + solntie_orig, gemtime_orig, 'b*', 'MarkerSize', 10); hold on; grid on;



