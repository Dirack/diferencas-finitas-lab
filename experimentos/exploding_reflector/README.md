# Modelo do refletor explosivo

> Este é um experimento de modelagem com diferenças finitas para exemplificar o modelo do refletor explosivo.

![Experimento numérico refletor explosivo](https://github.com/Dirack/diferencas-finitas-lab/blob/main/res/wav_com_borda.gif)

Neste experimento numérico, nós iremos gerar uma onda plana utilizando um conjunto de fontes pontuais (refletores explosivos) no centro do modelo.
A ideia é ilustrar o princípio de Huygens, da superposição de difrações para formar uma frente de onda plana, e o modelo do refletor explosivo.

Este experimento numético é reproduzido utilizando o core do pacote de processamento sísmico Madagascar 3.0. Estas dependências estão armazenaas nas
pastas bin, lib e include dentro do diretório 'experimentos', de modo a tornar este programa 'standalone'. Você pode compilar e rodar este experimento
no Linux Ubuntu utilizando 'make' ou rodar em um container do Docker, como descrito a seguir.

## Quick Start (Rodar o experimento no Linux com o Madagascar já instalado)

Você pode compilar este experimento no Linux Ubuntu como descrito a seguir. Porém, precisará do Madagascar 3.0 corretamente
instalado para rodar a modelagem com diferenças finitas e visualizar os resultados.
Para compilar o programa de modelagem 'Mfdmodeling.c' utilize o comando:

```
make fdmodeling.x
```

Para gerar o modelo de velocidades e rodar a modelagem utilize:

```
make experiment
```

Para vizualizar o experimento numérico utilize:

```
make view
```

E para rodar os 3 comandos acima em sequência utilize:

```
make
```

# Compilar e rodar o experimento em um container do docker

Primeiro, utilizando o Dockerfile deste diretório criamos a imagem docker com o Madagascar 3.0 instalado e configurado:

```
docker build -t fdexperimentos . -f Dockerfile
```

Depois de construir a imagem 'fdexperimentos', podemos criar um container do docker para fazer a compilação e rodar o experimento numérico
gerando os arquivos de saída utilizando o diretório atual como volume mapeando o diretório '/home/tryitondocker/experimentos' dentro do container:

```
docker run -ti -v$(pwd):/home/tryitondocker/experimentos fdexperimentos
```

Executamos o seguinte comando dentro do container e os arquivos de saída serão gerados:

```
make fdmodeling.x experiment
```

Daí é só sair do container e utilizar o comando 'make view' para visualizar o resultado do experimento numérico.

# Converter resultado de vpl para gif utilizando o utilitário do Madagascar

O pacote Madagascar possui o programa vpconvert para conversão de arquivos de imagem do Madagascar (.vpl) para diferentes formatos
de imagem. Se você tiver o pacote Madagascar instalado, poderá utilizar o utilitário para converter o resultado do experimento numérico,
armazenado em 'wav_com_borda.vpl', para uma imagem gif com o comando a seguir:

```
vpconvert format=gif wav_com_borda.vpl
```

Se desejar visualizar a partir do terminal do shell, utilize o comando a seguir:

```
xdg-open wav_com_borda.gif
```
