using Adchou_sample
using Documenter

DocMeta.setdocmeta!(Adchou_sample, :DocTestSetup, :(using Adchou_sample); recursive=true)

makedocs(;
    modules=[Adchou_sample],
    authors="romgoti <romain.angotti@sciencespo.fr> and contributors",
    sitename="Adchou_sample.jl",
    format=Documenter.HTML(;
        canonical="https://romgoti.github.io/Adchou_sample.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/romgoti/Adchou_sample.jl",
    devbranch="main",
)
