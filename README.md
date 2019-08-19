# RCC - Rocket Cooling Calculator

O RCC é a implementação numérica em Python de um modelo de refrigeração regenerativa, descrito por Foltran e Blavier (2018) para cálculo do perfil de temperaturas para motores foguete bi propelente. Também inclui uma funcionalidade de otimização da geometria dos canais através do método do gradiente descendente.
Desenvolvido em 2019 por Jefferson Bezerra na forma de um trabalho de conclusão de curso, do bacharelado em Engenharia Mecânica da Universidade Federal do Ceará, sob orientação do Prof. Claus Wehmann.

O programa está dispońivel na forma de uma biblioteca para Python (librcc) e uma interface gráfica (RCC GUI).
<img src="/figures/fig1.png" width="300"/>

A validação do código é realizada utilizando o motor foguete bi proplente L-75, desenvolvido pela Força Aérea Brasileira (FAB). Realizada através da comparação do fluxo de calor ao longo do motor com os resultados de Almeida e Shimote (1999) e de Foltran e Blavier (2018), além de um comparativo do perfil de temperaturas com o último trabalho citado.
<img src="/figures/fig2.png" width="300"/>

# Requisitos e recomendações
- Debian 9 (Stretch) x64 com XFCE
> Provavelmente funcionará em outras distribuições e em ambiente Windows.
- Python 3.5
> Não testado em outras versões.
- PySimpleGUI
> Necessário para utilizar a interface gráfica.
- NumPy
- SciPy
- pyCEA (modificado)
> Utilizar a versão modificada disponível no repositório.
- GFortran
- NASA Chemical Equilibrium with Applications (CEA)
> O CEA e o Fortran necessitam estar dentro da pasta "CEA+FORTRAN".

# Estrutura

|Arquivo       |Descrição                         |
|--------------|----------------------------------------------------|
|l75.csv       |Geometria do motor foguete L-75                     |
|l75.rcc       |Exemplo configuração salva para uso no rccgui.py    |
|librcc.py     |Biblioteca Rocket Cooling Calculator                |
|ogum.csv      |Geometria do motor foguete Ogum                     |
|ogum.rcc      |Exemplo configuração salva para uso no rccgui.py    |
|pyCEA.py      |Biblioteca pyCEA modificada para uso com o librcc.py|
|rccgui.py     |Interface gráfica para utilização do librcc.py      |

# Utilização

O programa pode ser utilizado na forma de biblioteca e integrado em outras aplicações ou através da interface gráfica.

# Referências

Para consulta completa das referências acesse o arquivo TCC.pdf (disponível em breve) na pasta documentos. O diretório também contém artigos utilizados nas referências do trabalho.