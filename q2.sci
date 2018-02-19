//[p,m,x_final,y_final] = ag(20,20,-5,-5,5,5,0.85,0.01,50);
//Chamada geral da função, recebe como parâmetro a quantidade de Indivíduos e de bits desejada, os limites em X e Y, as taxas de crossover e mutação, e a quantidade de gerações; retorna a população após todas as gerações, o valor mínimo de cada geração e os valores normalizados dos cromossomos para futuras avaliações dos resultados
function [populacaoFinal,minimo,x_final,y_final] = ag(qntIndividuos, bits, bottom_x, bottom_y, top_x, top_y, taxaCrossover, taxaMutacao, nGeracoes)
    printf("Iniciando Algoritmo Genético\n")
    printf("Intervalo de %d a %d em x\n",bottom_x,top_x)
    printf("Intervalo de %d a %d em y\n",bottom_y,top_y)
    printf("Taxa de Crossover = %.3f, Taxa de Mutação = %.3f\n",taxaCrossover,taxaMutacao)
    printf("População de %d Indivíduos com representação de %d bits\n",qntIndividuos,bits)
    printf("População Inicial:\n")
    //gera e exibe a população inicial
    Populacao = gerarPopulacaoInicial(qntIndividuos,bits,bottom_x,bottom_y,top_x,top_y)
    for j = 1:size(Populacao,1)
        x = normalizacao(Populacao(j,1),bottom_x,top_x)
        y = normalizacao(Populacao(j,2),bottom_y,top_y)
        printf("  %s , %s  =  (%f,%f)  |  f(x,y) = %f\n",Populacao(j,1),Populacao(j,2), x, y, func([x;y]))
    end
    //métodos aplicados a todas as gerações
    for i = 1:nGeracoes
        //gera a próxima população pelo torneio, cruzamento e mutação
        Populacao = torneio(Populacao,bottom_x,bottom_y,top_x,top_y)
        Populacao = crossover(Populacao,taxaCrossover,bottom_x,bottom_y,top_x,top_y)
        Populacao = mutacao(Populacao,taxaMutacao,bottom_x,bottom_y,top_x,top_y)
        //exibe a atual população
        printf("%da Geração:\n",i)
        for j = 1:size(Populacao,1)
            x = normalizacao(Populacao(j,1),bottom_x,top_x)
            y = normalizacao(Populacao(j,2),bottom_y,top_y)
            x_final(j,i) = x
            y_final(j,i) = y
            printf("  %s , %s = (%f,%f) | f(x,y) = %f | g(x,y) = %f\n",Populacao(j,1),Populacao(j,2), x,y, func([x;y]),(2*x^2 - 5*x*y + 2*y^2))
        end
        minimo($+1) = resultados(Populacao,bottom_x,bottom_y,top_x,top_y)
    end
    populacaoFinal = Populacao
endfunction

//Resultados: recebe como parâmetro a população de uma geração e retorna o menor valor da função presente na população
function minimo = resultados(populacao,bx,by,tx,ty)
    s1 = normalizacao(populacao(:,1),bottom_x,top_x)
    s2 = normalizacao(populacao(:,2),bottom_y,top_y)
    s = [s1,s2]
    f = func(s')
    minimo = min(f)
endfunction

function melhores = torneio(populacao,bottom_x,bottom_y,top_x,top_y)
    qnt = size(populacao,1)
    melhores = []
    for i = 1:qnt
        k = ceil(rand()*qnt)
        a_x = populacao(k,1)
        a_y = populacao(k,2)
        k = ceil(rand()*qnt)
        b_x = populacao(k,1)
        b_y = populacao(k,2)
        if valorFuncao(a_x,a_y,bottom_x,bottom_y,top_x,top_y) > valorFuncao(b_x,b_y,bottom_x,bottom_y,top_x,top_y) then
            melhores(i,1) = b_x
            melhores(i,2) = b_y
        else
            melhores(i,1) = a_x
            melhores(i,2) = a_y
        end
    end
endfunction

//Gerar Cromossomo: recebe como parâmetro o cromossomo X e de acordo com a restrição e os limites, retorna o cromossomo Y
function [y_bin,r] = gerarCromossomo(x_bin,bx,by,tx,ty)
    r = 0
    x = normalizacao(x_bin,bx,tx)
    //cria um polinômio que a partir da restrição uma vez que já temos o valor de x
    z = poly([2*x^2,-5*x,2],'y',"coeff")
    //as raízes desse polinômio são os possíveis valores para y
    w = roots(z)
    a = evstr(string(w(1)))
    b = evstr(string(w(2)))
    //escolhe qualquer uma das raízes contanto que esteja dentro dos limites e não seja complexo
    if (a > by & a < ty & imag(a) == 0)
        if (b > by & b < ty & imag(b) == 0) then
            y = w(ceil(rand()*2))
        else
            y = w(1)
        end
    else
        if (b > by & b < ty & imag(b) == 0) then
            y = w(2)
        else
            r = 1
        end
    end
    s = evstr(string(y))
    //retorna a string binária
    y_bin = normalizacao_inv(s,by,ty,length(x_bin))
endfunction

//Crossover: recebe como parâmetro a população ordenada de acordo com a seleção e retorna os filhos na mesma ordem mas após cruzamento a cada dois indivíduos
function filhos = crossover(melhores,taxaCross,bx,by,tx,ty)
    qnt = size(melhores,1)
    filhos = melhores
    //avança a cada dois indivíduos
    for i = 1:2:qnt
        //testa aleatoriamente se há cruzamento no cromossomo X
        if (rand() <= taxaCross) then
            repeat = 1
            //repete até que tenha X e Y de acordo com a restrição
            while (repeat)
                //define a região de corte e separa a string binária dos cromossomos X dos dois indivíduos
                corte = ceil(rand()*(length(filhos(i,1))-1))
                str1 = strsplit(filhos(i,1),corte)
                str2 = strsplit(filhos(i+1,1),corte)
                //faz a concatenação das cabeças com as caudas dos pais
                p = str1(1)+str2(2)
                q = str2(1)+str1(2)
                //gera o cromossomo Y do primeiro filho
                [p2,repeat] = gerarCromossomo(p,bx,by,tx,ty)
                if (~repeat) then
                    //e do segundo
                    [q2,repeat] = gerarCromossomo(q,bx,by,tx,ty)
                end
                //se não há um Y adequado para a restrição dado o X, repete o cruzamento
            end
            //encontrou um ponto que atende à restrição
            filhos(i,1) = p
            filhos(i,2) = p2
            filhos(i+1,1) = q
            filhos(i+1,2) = q2
        end
    end
endfunction

//Mutação: recebe como parâmetro os filhos que foram gerados no cruzamento e retorna a população depois das mutações serem aplicadas
function novos = mutacao(filhos,taxaMut,bx,by,tx,ty)
    novos = filhos
    //percorre toda a população, esceto o último indivíduo
    for i = 1:size(filhos,1)-1
        //percorre todos os bits(genes) do crmossomo X do indivíduo
        for j = 1:length(filhos(i,1))
            repeat = 1
            //repete até que tenha X e Y de acordo com a restrição
            while(repeat)
                repeat = 0;
                //testa se haverá mutação para esse bit
                if (rand() <= taxaMut) then
                    //separa o cromossomo, altera o bit e concatena novamente
                    a = part(novos(i,1),[1:(j-1)])
                    b = part(novos(i,1),j)
                    c = part(novos(i,1),[(j+1):length(novos(i,1))])
                    str = [a,b,c]
                    if (str(2) == '0') then str(2) = '1'
                    else str(2) = '0'
                    end
                    p = str(1)+str(2)+str(3)
                    //gera o cromossomo Y
                    [p2,repeat] = gerarCromossomo(p,bx,by,tx,ty)
                    if (~repeat) then
                        novos(i,1) = p
                        novos(i,2) = p2
                    end
                end
            end
        end
    end
endfunction

function apt = valorFuncao(ind_x,ind_y,bottom_x,bottom_y,top_x,top_y)
    x = normalizacao(ind_x,bottom_x,top_x)
    y = normalizacao(ind_y,bottom_y,top_y)
    apt = func([x;y])
endfunction

//Gerar População Inicial: gera aleatoriamente uma populacao inicial de acordo com a quantidade de indivíduos e de bits desejadas que são passadas como parâmetros
function Inicial = gerarPopulacaoInicial(qntInd,qntBits,bx,by,tx,ty)
    for i = 1:qntInd
        repeat = 1
        //repete até que tenha X e Y de acordo com a restrição
        while (repeat)
            //gera todos os bits do cromossomo X do indivíduo em questão
            p = string(int(rand(1, qntBits)*2))
            x_bin = strcat(p(1,:))
            //gera o cromossomo Y
            [s,repeat] = gerarCromossomo(x_bin,bx,by,tx,ty)
        end
        Inicial(i,1) = x_bin
        Inicial(i,2) = s
    end
endfunction

function dec = normalizacao(bin, b, t)
    dec = b+(t-b) .* bin2dec(bin) ./ ((2^length(bin)) -1)
endfunction

function bin = normalizacao_inv(dec, b, t, n_bits)
    bin = dec2bin((dec-b) .* ((2^n_bits)-1) ./ (t-b), n_bits)
endfunction

function f = func(p)
    x = p(1,:)
    y = p(2,:)
    f = (x-1).^2 + (y-1).^2
    f = f'
endfunction