function [G1S,G2M] = get_ccscore(X,genes,do_plot)
% returns G1S and G2M score; X is expression matrix ([patient,genes]),
%   genes is list of symbols.
data=X;
cc_file = 'cc_genes.csv'; % based on the list described by Tirosh et al.

if ~strcmp(cc_file, '');
    
    M = readtable(cc_file,'ReadVariableNames',1);

    G1S_genes = table2cell(M(:,1))';
    G2M_genes = table2cell(M(:,2))';

    G1S_genes = G1S_genes(~strcmp(G1S_genes,''));
    G2M_genes = G2M_genes(~strcmp(G2M_genes,''));

    G1Si=[];
    G2Mi=[];

    disp('-- Indexing G1S');
    G1S_not_found = 0;
    for i = 1:numel(G1S_genes);
        index = find(strcmpi(genes,G1S_genes{i}));
        if ~ isempty(index);
            if length(index) == 1;
                G1Si(end+1) = index;
            else
                disp([G1S_genes{i} ' is found more than once']);
            end
        else
            disp([G1S_genes{i} ' cannot be found']);
            G1S_not_found = G1S_not_found +1;
        end
    end

    disp('-- Indexing G2M');
    G2M_not_found = 0;
    for i = 1:numel(G2M_genes);
        index = find(strcmpi(genes,G2M_genes{i}));
        if ~ isempty(index);
            if length(index) == 1;
                G2Mi(end+1) = index;
            else
                disp([G2M_genes{i} ' is found more than once']);
            end
        else
            disp([G2M_genes{i} ' cannot be found']);
            G2M_not_found = G2M_not_found +1;
        end
    end
    disp(newline)
    fprintf('%.2f%% G1S genes used\n',(numel(G1S_genes)-G1S_not_found)/numel(G1S_genes)*100);
    fprintf('%.2f%% G2M genes used\n',(numel(G2M_genes)-G2M_not_found)/numel(G2M_genes)*100);

    % phase score
    for i = 1:size(data,1);
        G1S(i) = mean(data(i,G1Si));
        G2M(i) = mean(data(i,G2Mi));
    end

    G1S=zscore(G1S);
    G2M=zscore(G2M);
    
    if do_plot
        % Plot
        figure()
        scatter(G1S,G2M,5,'fill')
        xlabel('G1/S score')
        ylabel('G2/M score')
        title('cell cycle phase')
    end
end


end

