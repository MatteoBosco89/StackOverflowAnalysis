%% Studio StackOverflow
% carichiamo anno di interesse
lista = fopen("filtered-2015.txt", "rt");
C = textscan(lista, '%d %d %d');
fclose(lista);
s = C(1);
d = C(2);
t = C(3);
n = size(s{1});
n = n(1);
% inizializzo sorgente e destinazione
sor = zeros(1, n);
des = zeros(1, n);
% mappiamo il grafo
assMap = containers.Map();
% carichiamo assMap con valori random
% la chiave è importante, ad ogni chiave è associato un utente reale
for i=1:n
    assMap(""+s{1}(i)) = i;
    assMap(""+d{1}(i)) = i;
end
% costruzione etichettamento nodi
% prendo le chiavi di assMap
k = keys(assMap);
% con i indico il nodo su matlab
% con k{i} indico l'utente (etichetta)
% gli utenti in k non sono ordinati
for i = 1:length(assMap)
    assMap(k{i}) = i;
end
% utilizzo mappa come wrapper
% n esprime il numero di archi
for i = 1:n
    sor(1, i) = assMap(""+s{1}(i));
    des(1, i) = assMap(""+d{1}(i));
end
% quindi, per avere l'utente a partire dal numero di nodo di MatLab si può
% semplicemente utilizzare utente = k{nodo_matlab}

% assegno ai nodi di MatLab il numero dell'utente del dataset
% assegnando il nome al nodo
nodenames = zeros(1, length(assMap));
for i=1:length(assMap)
    nodenames(i) = ""+k{i}; 
end

% plot grafo con tutte le componenti con e senza self loops
Goriginale = digraph(sor, des, "omitselfloops");

%% Eliminazione archi duplicati e calcolo pesi archi
Goriginale.Edges.Weigth = ones(height(Goriginale.Edges), 1);

G = simplify(Goriginale, 'sum', 'AggregationVariables', {'Weigth'});

%% Assegnazione nomi utenti a nodi
G.Nodes.Name = k';
% vettore dei nomi k'
Gl = digraph(sor, des);

%% costruisco matrice di adiacenza
A = adjacency(G);
% A = full(adjacency(G));
% plot con etichette ai nodi
% figure, plot(G,'NodeLabel',nodenames);
% plot normali
% figure, plot(G);
% figure, plot(Gl);

%% Analisi Temporale grado
indeg = indegree(Goriginale);
outdeg = outdegree(Goriginale);
deg = indeg + outdeg;

% analizziamo evoluzione durante gli anni
meandeg = mean(deg);

%% Componente Connessa
% rimuoviamo gli utenti con meno di 50 interazioni
degC = deg <= 49;
CG = rmnode(G, find(degC == 1));
ACG = adjacency(CG);
% ripulisco le componenti
[bin, binsize] = conncomp(CG, 'Type', 'strong');
% componente connessa più grande
comp1 = max(binsize);
idx = binsize(bin) == comp1;
% seconda componente connessa più grande
[aa, aaa] = max(binsize);
binsize(aaa) = 0;
comp2 = max(binsize);
% prendiamo in considerazione la più grande trovata
CG = subgraph(CG, idx);
Aconn = adjacency(CG);
Afilt = adjacency(CG);


%% Creazione matrice adiacenza pesata
Apesata = adjacency(CG, 'weighted');

%% Grafo finale
Gcomp = CG;
% numero nodi grafo
n = height(Gcomp.Nodes);

%% Grado e Closeness
D = distances(Gcomp);
% verifica fortemente connesso 
Dv = D + eye(n);
if(min(min(Dv)) == 1)
    disp 'Fortemente Connesso'
end


Dinv = 1./D;
Dinv(isinf(Dinv)) = 0;

% calcoliamo le interazioni
gradoin = indegree(Gcomp);
gradoout = outdegree(Gcomp);
grado = gradoin + gradoout;
%figure, bar(gradoin/max(gradoin));
%figure, bar(gradoout/max(gradoout));
% determiniamo i nodi al centro della rete e ai margini
cin = sum(Dinv, 1)/(n-1);
cout = sum(Dinv, 2)/(n-1);
%figure, 
%h = bar(cin);
%hold on;
%h = bar(cout');
%legend('In-Closeness','Out-Closeness');


%% STUDIO CONNESSIONE DEL GRAFO
% eliminiamo i nodi dal grafo che hanno grado più alto e vediamo se il
% grafo rimane ancora fortemente connesso
Gverifica = Gcomp;
utenti_importanti_nodo = [];
utenti_importanti_reale = [];
[bin, binsize] = conncomp(Gverifica);
while(max(binsize) == height(Gverifica.Nodes))
    gradoinV = indegree(Gverifica);
    gradooutV = outdegree(Gverifica);
    gradoV = gradoinV + gradooutV;
    [gmax, utente] = max(gradoV);
    utenti_importanti_nodo = [utenti_importanti_nodo; utente];
    utenti_importanti_reale = [utenti_importanti_reale; Gverifica.Nodes(utente, 1)];
    Gverifica = rmnode(Gverifica, utente);
    [bin, binsize] = conncomp(Gverifica);
end


%% Densità del grafo m/n*n-1
% grafo di partenza 
nG = height(G.Nodes);
mG = height(G.Edges);
DensityG = mG / (nG*nG - 1);
% componente fortemente connessa
nCG = height(CG.Nodes);
mCG = height(CG.Edges);
DensityCG = mCG / (nCG*nCG - 1);

%% Verifica importanza strategica nodo
% calcoliamo se il nodo con massima betweenness è proprio quello eliminato
% prima
betweenness = centrality(Gcomp, 'betweenness');
figure, bar(betweenness/max(betweenness));
[bmax, bpos] = max(betweenness);

%% grado in ingresso e uscita distribuzione
gin = zeros(length(gradoin), 1);
for i=1:length(gradoin)
    gin(i) = sum(gradoin == i)/length(gradoin);
end
gout = zeros(length(gradoout), 1);
for i=1:length(gradoout)
    gout(i) = sum(gradoout == i)/length(gradoout);
end
figure, bar([gin gout]);
legend('Indegree', 'Outdegree');

%% calcolo hub e autority 
% gli hub sono gli utenti che hanno risposto alle domande più importanti
% le autority sono gli utenti che hanno posto le domande che hanno ricevuto
% più risposte e possono essere viste come domande con argomenti più
% importanti
% controlliamo se gli hub sono gli stessi nodi scollegati prima
% possiamo vedere se davvero il grado in uscita è più importante del grado
% in ingresso
hub_ranks = centrality(Gcomp, 'hubs');
auth_ranks = centrality(Gcomp, 'authorities');
Gcomp.Nodes.Hubs = hub_ranks;
Gcomp.Nodes.Authorities = auth_ranks;
% verifichiamo se l'hub maggiore è l'utente trovato prima
[hmax, hpos] = max(hub_ranks);
% calcoliamo l'authority maggiore
[amax, apos] = max(auth_ranks);

%% CODICE BOXPLOT GRADI PER ANNO
x = [deg2009; deg2010; deg2011; deg2012; deg2013; deg2014; deg2015];
g = [zeros(length(deg2009), 1); ones(length(deg2010), 1); 2*ones(length(deg2011), 1); 3*ones(length(deg2012), 1); 4*ones(length(deg2013), 1); 5*ones(length(deg2014), 1); 6*ones(length(deg2015), 1)];
mediana = [median(deg2009) median(deg2010) median(deg2011) median(deg2012) median(deg2013) median(deg2014) median(deg2015)];
quartile10 = [quantile(deg2009, 0.10) quantile(deg2010, 0.10) quantile(deg2011, 0.10) quantile(deg2012, 0.10) quantile(deg2013, 0.10) quantile(deg2014, 0.10) quantile(deg2015, 0.10)];
quartile25 = [quantile(deg2009, 0.25) quantile(deg2010, 0.25) quantile(deg2011, 0.25) quantile(deg2012, 0.25) quantile(deg2013, 0.25) quantile(deg2014, 0.25) quantile(deg2015, 0.25)];
quartile75 = [quantile(deg2009, 0.75) quantile(deg2010, 0.75) quantile(deg2011, 0.75) quantile(deg2012, 0.75) quantile(deg2013, 0.75) quantile(deg2014, 0.75) quantile(deg2015, 0.75)];
quartile90 = [quantile(deg2009, 0.90) quantile(deg2010, 0.90) quantile(deg2011, 0.90) quantile(deg2012, 0.90) quantile(deg2013, 0.90) quantile(deg2014, 0.90) quantile(deg2015, 0.90)];
anni = [2009 2010 2011 2012 2013 2014 2015];
boxplot(x, g,'Labels',{'2009','2010', '2011','2012', '2013','2014', '2015'});
xlabel('Anni');
ylabel('Grado');

figure, plot(anni, mediana, 'r');
hold on
plot(anni, quartile10, 'black');
plot(anni, quartile25, 'g');
plot(anni, quartile75, 'y');
plot(anni, quartile90);
legend('Median', '10th percentile', '25th percentile', '75th percentile', '90th percentile');

%% TROVARE NODO IN BASE AL NOME UTENTE
% utenti 22656, 548225
nomeUtente = "";
nodo_u = find(ismember(CG.Nodes.Name, nomeUtente) == 1);




