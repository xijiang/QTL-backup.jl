"""
    function gs_accuracy(Nₑ, L, k, h², T, n; ϵ = 1e-4)
Given a selective population's effective size, `Nₑ`,
length of genome, `L`, in Morgen,
number of chromosome, `k`,
trait heritability, `h²`,
number of phenotyped individuals, `T`,
number of available markers,
this function returns the prediction accuracy with genomic selection.

Here, accuracy is the correlation between the true breeding values
and the GEBV.
- Refer doi: 10.1146/annurev-animal-031412-103-705
"""
function gs_accuracy(Nₑ, L, k, h², T, n; ϵ = 1e-4)
    Mₑ = 2Nₑ * L * k / log(2Nₑ)
    #β = σᵤ² / (σᵤ² + σₐ²)
    β = n / (n + Mₑ)
    θ = T * h² * β / Mₑ
    @show θ
    rp, r² = 0, .5
    while abs(r² - rp) > ϵ
        rp = r²
        r² = θ * β / (θ + 1 - h² * r²)
        println(r²)
    end
    r²
end
