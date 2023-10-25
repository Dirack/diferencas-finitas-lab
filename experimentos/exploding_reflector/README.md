# Modelo do refletor explosivo

Este é um experimento de modelagem com diferenças finitas para exemplificar o modelo do refletor explosivo.
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

Primeiro, utilizando o Dockerfile deste diretório criamos a imagem docker com o Madagascar instalado e configurado:

```
docker build -t fdexperimentos . -f Dockerfile
```

Depois de construir a imagem, podemos criar um container do docker para fazer a compilação e rodar o experimento numérico
gerando os arquivos de saída:

```
docker run -ti -v$(pwd):/home/tryitondocker/experimentos fdexperimentos
```

Executamos o seguinte comando dentro do container e os arquivos de saída serão gerados:

```
make fdmodeling.x experiment
```

Daí é só sair do container e utilizar o comando 'make view' para visualizar o resultado do experimento numérico.
