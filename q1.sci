//[p,m,x_final,y_final] = ag(75,20,-500,-500,500,500,0.85,0.005,50)
//Chamada geral da função, recebe como parâmetro a quantidade de Indivíduos e de bits desejada, os limites em X e Y, as taxas de crossover e mutação, e a quantidade de gerações; retorna a população após todas as gerações, o valor mínimo de cada geração e os valores normalizados dos cromossomos para futuras avaliações dos resultados
function [populacaoFinal,minimo,x_final,y_final] = ag(qntIndividuos, bits, bottom_x, bottom_y, top_x, top_y, taxaCrossover, taxaMutacao, nGeracoes)
    printf("Iniciando Algoritmo Genético\n")
    printf("Intervalo de %d a %d em x\n",bottom_x,top_x)
    printf("Intervalo de %d a %d em y\n",bottom_y,top_y)
    printf("Taxa de Crossover = %.3f, Taxa de Mutação = %.3f\n",taxaCrossover,taxaMutacao)
    printf("População de %d Indivíduos com representação de %d bits\n",qntIndividuos,bits)
    printf("População Inicial:\n")
    //gera e exibe a população inicial
    Populacao = gerarPopulacaoInicial(qntIndividuos,bits)
    for j = 1:size(Populacao,1)
        x = normalizacao(Populacao(j,1),bottom_x,top_x)
        y = normalizacao(Populacao(j,2),bottom_y,top_y)
        printf("  %s , %s  =  (%f,%f)  |  f(x,y) = %f\n",Populacao(j,1),Populacao(j,2), x, y, func([x;y]))
    end
    //métodos aplicados a todas as gerações
    for i = 1:nGeracoes
        //aloca o indivíduo otimizado por Newton para colocálo na próxima geração
        Melhor = Newton(Populacao,bottom_x,bottom_y,top_x,top_y)
        //gera a próxima população pelo torneio, cruzamento e mutação
        Populacao = torneio(Populacao,bottom_x,bottom_y,top_x,top_y)
        Populacao = crossover(Populacao,taxaCrossover)
        Populacao = mutacao(Populacao,taxaMutacao)
        //coloca o indivíduo otimizado por Newto na nova população
        Populacao(qntIndividuos,1) = Melhor(1)
        Populacao(qntIndividuos,2) = Melhor(2)
        //exibe a atual população
        printf("%da Geração:\n",i)
        for j = 1:size(Populacao,1)
            x = normalizacao(Populacao(j,1),bottom_x,top_x)
            y = normalizacao(Populacao(j,2),bottom_y,top_y)
            x_final(j,i) = x
            y_final(j,i) = y
            printf("  %s , %s  =  (%f,%f)  |  f(x,y) = %f\n",Populacao(j,1),Populacao(j,2), x,y, func([x;y]))
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

//Newton: recebe toda a população como parâmetro e retorna um indivíduo (em binário) resultado do método de Newton aplicado ao indivíduo mais apto da população
function melhor = Newton(populacao,bottom_x,bottom_y,top_x,top_y)
    qnt = size(populacao,1)
    //normaliza todos os indivíduos, calcula seus valores na função e acha o menor
    s1 = normalizacao(populacao(:,1),bottom_x,top_x)
    s2 = normalizacao(populacao(:,2),bottom_y,top_y)
    s = [s1,s2]
    f = func(s')
    [minimo,index] = min(f)
    //copia o indivíduo mais apto, normaliza-o e acha seu valor na funcao, alocando também em variáveis auxiliares
    bin_x_min = [populacao(index,1);populacao(index,2)]
    n_bits = length(bin_x_min(1))
    p = normalizacao(bin_x_min(1), bottom_x, top_x)
    q = normalizacao(bin_x_min(2), bottom_y, top_y)
    x_min = [p;q]
    f_k = func(x_min)
    x_min2 = x_min
    f_k2 = f_k
    repeat = 1
    while (repeat)
        f_k = f_k2
        x_min = x_min2
        //encontra a Hessiana e o Gradiente da função no ponto em questão e calcula o próximo ponto
        [g,H] = derivative(func,x_min,H_form='blockmat')
        g = g'
        x_min2 = x_min - inv(H)*g
        f_k2 = func(x_min)
        //testa se os pontos estão fora dos limites ou se a função não está mais sendo minimizada, em caso positivo o ponto anterior é o ótimo
        if ((x_min2(1) < bottom_x) | (x_min2(1) > top_x)) then
            repeat = 0
        end
        if ((x_min2(2) < bottom_y) | (x_min2(2) > top_y)) then
            repeat = 0
        end
        if (f_k2 > f_k) then
            repeat = 0
        end
    end
    //retorna o ótimo em binário
    melhor(1) = normalizacao_inv(x_min(1),bottom_x,top_y,n_bits)
    melhor(2) = normalizacao_inv(x_min(2),bottom_x,top_y,n_bits)
endfunction

//Torneio: recebe como parâmetro toda a população e retorna uma população modificada com os indivíduos na ordem em que foram selecionados
function melhores = torneio(populacao,bottom_x,bottom_y,top_x,top_y)
    qnt = size(populacao,1)
    melhores = []
    //a quantidade de indivíduos selecionados para reprodução é a quantidade de indivíduos na população menos um indivíduo que será escolhido para aplicar Newton
    for i = 1:qnt-1
        //escolhe aleatoriamente dois individuos em toda a populacao
        k = ceil(rand()*qnt)
        a_x = populacao(k,1)
        a_y = populacao(k,2)
        k = ceil(rand()*qnt)
        b_x = populacao(k,1)
        b_y = populacao(k,2)
        //escolhe entre os dois de acordo com o valor deles na função
        a = valorFuncao(a_x,a_y,bottom_x,bottom_y,top_x,top_y)
        b = valorFuncao(b_x,b_y,bottom_x,bottom_y,top_x,top_y)
        if a > b then
            melhores(i,1) = b_x
            melhores(i,2) = b_y
        else
            melhores(i,1) = a_x
            melhores(i,2) = a_y
        end
    end
endfunction

//Crossover: recebe como parâmetro a população ordenada de acordo com a seleção e retorna os filhos na mesma ordem mas após cruzamento a cada dois indivíduos
function filhos = crossover(melhores,taxaCross)
    qnt = size(melhores,1)
    filhos = melhores
    //avança a cada dois indivíduos mas não faz cruzamento no último que será o resultado do método de Newton
    for i = 1:2:qnt-1
        //testa aleatoriamente se há cruzamento no cromossomo X
        if (rand() <= taxaCross) then
            //define a região de corte e separa a string binária dos cromossomos X dos dois indivíduos
            corte = ceil(rand()*(length(filhos(i,1))-1))
            str1 = strsplit(filhos(i,1),corte)
            str2 = strsplit(filhos(i+1,1),corte)
            //faz a concatenação das cabeças com as caudas dos pais
            filhos(i,1) = str1(1)+str2(2)
            filhos(i+1,1) = str2(1)+str1(2)
        end
        
        //faz o mesmo para o cromossomo Y
        if (rand() <= taxaCross) then
            corte = ceil(rand()*(length(filhos(i,1))-1))
            str1 = strsplit(filhos(i,2),corte)
            str2 = strsplit(filhos(i+1,2),corte)
            filhos(i,2) = str1(1)+str2(2)
            filhos(i+1,2) = str2(1)+str1(2)
        end
    end
endfunction

//Mutação: recebe como parâmetro os filhos que foram gerados no cruzamento e retorna a população depois das mutações serem aplicadas
function novos = mutacao(filhos,taxaMut)
    novos = filhos
    //percorre toda a população, esceto o último indivíduo
    for i = 1:size(filhos,1)-1
        //percorre todos os bits(genes) do crmossomo X do indivíduo
        for j = 1:length(filhos(i,1))
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
                novos(i,1) = str(1)+str(2)+str(3)
            end
        end
        
        //faz o mesmo para o cromossomo Y
        for j = 1:length(filhos(i,2))
            if (rand() <= taxaMut) then
                a = part(novos(i,1),[1:(j-1)])
                b = part(novos(i,1),j)
                c = part(novos(i,1),[(j+1):length(novos(i,1))])
                str = [a,b,c]
                if (str(2) == '0') then str(2) = '1'
                else str(2) = '0'
                end
                novos(i,2) = str(1)+str(2)+str(3)
            end
        end
    end
endfunction

//Valor na Função: recebe como parâmetro os cromossomos X e Y de um indivíduo e retorna o valor dele na função
function apt = valorFuncao(ind_x,ind_y,bottom_x,bottom_y,top_x,top_y)
    x = normalizacao(ind_x,bottom_x,top_x)
    y = normalizacao(ind_y,bottom_y,top_y)
    apt = func([x;y])
endfunction

//Gerar População Inicial: gera aleatoriamente uma populacao inicial de acordo com a quantidade de indivíduos e de bits desejadas que são passadas como parâmetros
function Inicial = gerarPopulacaoInicial(qntInd,qntBits)
    //gera todos os bits
    p = string(int(rand(qntInd, qntBits)*2))
    q = string(int(rand(qntInd, qntBits)*2))
    //concatena os bits de um cromossomo
    for i = 1:qntInd
        Inicial(i,1) = strcat(p(i,:))
        Inicial(i,2) = strcat(q(i,:))
    end
endfunction

//Normalização: recebe como parâmetro os limites e um cromossomo (ou matriz de cromossomos) e retorna o seu numero em decimal correspondente ao intervalo especificado
function dec = normalizacao(bin, b, t)
    dec = b+(t-b) .* bin2dec(bin) ./ ((2^length(bin)) -1)
endfunction

//Normalização Inversa: transforma o número decimal dentro do intervalo especificado no string binário do cromossomo
function bin = normalizacao_inv(dec, b, t, n_bits)
    bin = dec2bin((dec-b) .* ((2^n_bits)-1) ./ (t-b), n_bits)
endfunction

//Função a ser minimizada
function w11 = func(p)
    x = p(1,:)
    y = p(2,:)
    z=-x.*sin(sqrt(abs(x)))-y.*sin(sqrt(abs(y)));
    x=x/250;
    y=y/250;
    r=100*(y-x.^2).^2+(1-x).^2;
    w3=r-z;
    w4=sqrt(r.^2+z.^2);
    w11=w3.^2+w4;
    //w11 = -w11
    w11 = w11'
endfunction