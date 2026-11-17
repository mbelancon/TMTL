"""Translations of the user-facing strings produced by the Python side.

The code is written in English. Messages that reach the interface (model titles,
descriptions, parameter labels, θ definitions, fit log, errors) also have a
Portuguese version, kept here in the dictionary ``PT`` keyed by the English text
(or by the English *template*, for messages with fields such as ``{chi:.3e}``).

* ``Msg(template, **fields)`` is a ``str`` with the English text that remembers
  its template and fields, so it can be translated later with ``tr``.
* ``tr(text, lang)`` returns the text in ``lang`` ("en" or "pt"); unknown texts
  are returned unchanged (in English).

The interface's own strings (buttons, labels, plot axes) live in
``tlfit/static/index.html``.
"""
from __future__ import annotations

LANGS = ("en", "pt")


class Msg(str):
    """English message that keeps its template and fields for translation."""

    def __new__(cls, template: str, **fields):
        s = super().__new__(cls, template.format(**fields) if fields else template)
        s.template, s.fields = template, fields
        return s


def norm_lang(lang) -> str:
    return "pt" if str(lang or "").lower().startswith("pt") else "en"


def tr(text, lang="en") -> str:
    """Translates an English text (or a Msg) into ``lang``."""
    if norm_lang(lang) == "en":
        return str(text)
    if isinstance(text, Msg):
        tpl = PT.get(text.template)
        if tpl is None:
            return str(text)
        fields = {k: tr(v, lang) if isinstance(v, Msg) else v for k, v in text.fields.items()}
        return tpl.format(**fields)
    return PT.get(str(text), str(text))


def both(text) -> dict:
    """{"en": ..., "pt": ...} for a string shown in the interface."""
    return {lang: tr(text, lang) for lang in LANGS}


PT: dict[str, str] = {
    # ------------------------------------------------------------------ parameter labels
    "D sample": "D amostra",
    "θ sample": "θ amostra",
    "θ fluid": "θ fluido",
    "L thickness": "L espessura",
    "l thickness": "l espessura",
    "D fluid": "D fluido",
    "c sample": "c amostra",
    "c fluid": "c fluido",
    "ρ sample": "ρ amostra",
    "ρ fluid": "ρ fluido",
    "n sample": "n amostra",
    "n fluid": "n fluido",
    "α_T expansion": "α_T expansão",
    "dn/dT sample": "dn/dT amostra",
    "I0 (scale)": "I0 (escala)",
    "τ (lifetime)": "τ (vida)",
    "h (excitation rate)": "h (taxa de excitação)",
    "A_e absorption": "A_e absorção",
    "t_off (0 = no switch-off)": "t_off (0 = sem desligar)",
    "β expansion": "β expansão",
    "dn/dT air": "dn/dT ar",
    "D air": "D ar",
    "dimensionless": "adimensional",
    # ------------------------------------------------------------------ presets
    "Air": "Ar",
    "Water": "Água",
    # ------------------------------------------------------------------ model titles
    "Shen TL (classical)": "TL de Shen (clássico)",
    "TL with photobleaching (Herculano 2013)": "TL com fotobranqueamento (Herculano 2013)",
    "TL + population lens (Silva 2013)": "TL + lente de população (Silva 2013)",
    "TM LAM (low absorption)": "TM LAM (baixa absorção)",
    "TM LAM, finite sample": "TM LAM, amostra finita",
    "TM BAM (Beer's law)": "TM BAM (lei de Beer)",
    "TM+TL-Fluid, HAM on/off (approximation: no feedback from the fluid)":
        "TM+TL-Fluid, HAM on/off (aproximação: sem retorno do fluido)",
    "TM+TL-Fluid, HAM (full coupling, Flizikowski 2020)":
        "TM+TL-Fluid, HAM (acoplamento completo, Flizikowski 2020)",
    # ------------------------------------------------------------------ model descriptions
    "Thermal lens with sample–fluid heat coupling (Malacarne et al., Appl. Spectrosc. 65, 99, 2011). "
    "Equivalent to tllam.c, with the fluid term corrected.":
        "Lente térmica com acoplamento térmico amostra–fluido (Malacarne et al., Appl. Spectrosc. 65, 99, "
        "2011). Equivalente a tllam.c, com o termo do fluido corrigido.",
    "Closed-form expression of Shen's aberrant model for the thermal lens in the mode-mismatched "
    "configuration (thin sample, low absorption, small θ). t_c = w0e²/(4D). Baesso, Shen & Snook, "
    "J. Appl. Phys. 75, 3732 (1994), eq. 10.":
        "Expressão fechada do modelo aberrante de Shen para lente térmica no modo descasado (amostra fina, "
        "baixa absorção, θ pequeno). t_c = w0e²/(4D). Baesso, Shen & Snook, J. Appl. Phys. 75, 3732 "
        "(1994), eq. 10.",
    "Thermal lens of a liquid whose absorption decays by a photochemical reaction (pseudo-first order, "
    "rate K_T) down to ε̄·β. Fits θ, K_T and ε̄; fix D with a measurement of the pure solvent (K_T = 0). "
    "Herculano et al., J. Phys. Chem. B 117, 1932 (2013), eqs. 3–6.":
        "Lente térmica de líquido com absorção que decai por reação fotoquímica (pseudo-primeira ordem, "
        "taxa K_T) até ε̄·β. Ajusta θ, K_T e ε̄; fixe D com uma medida do solvente puro (K_T = 0). "
        "Herculano et al., J. Phys. Chem. B 117, 1932 (2013), eqs. 3–6.",
    "Thermal lens plus excited-state population lens (lifetime τ, excitation rate h = I0(σ_abs+σ_se)/E "
    "at the centre, absorption A_e, thickness l). Silva et al., Appl. Phys. Lett. 102, 141910 (2013), "
    "eqs. 1–3. Fits θ_th, θ_PL, τ and D; m and V come from the setup, h, A_e and l are known.":
        "Lente térmica mais lente de população do estado excitado (vida τ, taxa de excitação "
        "h = I0(σ_abs+σ_se)/E no centro, absorção A_e, espessura l). Silva et al., Appl. Phys. Lett. 102, "
        "141910 (2013), eqs. 1–3. Ajusta θ_th, θ_PL, τ e D; m e V são da montagem, h, A_e e l são "
        "conhecidos.",
    "Thermal mirror for a semi-infinite sample with a source uniform in z (A_e·thickness ≪ 1). "
    "LAM = Low Absorption Model. Sato et al., J. Appl. Phys. 104, 053520 (2008), eq. 32. "
    "θ_TM = −P_e A_e α_T(1+ν)φ/(λ_p k), in 1/m.":
        "Espelho térmico para amostra semi-infinita com fonte uniforme em z (A_e·espessura ≪ 1). "
        "LAM = Low Absorption Model. Sato et al., J. Appl. Phys. 104, 053520 (2008), eq. 32. "
        "θ_TM = −P_e A_e α_T(1+ν)φ/(λ_p k), em 1/m.",
    "Low-absorption thermal mirror in a sample of thickness l, with both faces free. Astrath et al., "
    "J. Opt. Soc. Am. B 28, 1735 (2011), eq. 10. Tends to the semi-infinite LAM for large l.":
        "Espelho térmico de baixa absorção em amostra de espessura l, com as duas faces livres. Astrath "
        "et al., J. Opt. Soc. Am. B 28, 1735 (2011), eq. 10. Tende ao LAM semi-infinito para l grande.",
    "Thermal mirror with absorption following Beer's law, valid for any A_e (semi-infinite sample); "
    "BAM = Beer Absorption Model. Sato et al., J. Appl. Phys. 104, 053520 (2008), eq. 33–34. Tends to "
    "LAM (A_e → 0) and to HAM (A_e → ∞, with θ_HAM = θ_TM/A_e).":
        "Espelho térmico com absorção pela lei de Beer, válido para qualquer A_e (amostra semi-infinita); "
        "BAM = Beer Absorption Model. Sato et al., J. Appl. Phys. 104, 053520 (2008), eq. 33–34. Tende ao "
        "LAM (A_e → 0) e ao HAM (A_e → ∞, com θ_HAM = θ_TM/A_e).",
    "HAM thermal mirror (High Absorption Model) with the excitation switched off at t = t_off, plus the "
    "thermal lens of the fluid (air), active if dn/dT ≠ 0. Equivalent to tmhamtl-onoff.c: approximation "
    "in which the sample surface is free; the fluid receives heat and forms the TL, but does not cool "
    "the sample back. For the full coupling (Flizikowski et al., Opt. Express 28, 7116, 2020) use "
    "TM+TL-Fluid, HAM (full coupling): in the three metals of the paper the difference is 1e-4 to 1e-3 "
    "of the signal with the laser on and reaches ~1% right after switch-off (stainless steel in air).":
        "Espelho térmico HAM (High Absorption Model, alta absorção) com excitação desligada em t = t_off, "
        "mais a lente térmica do fluido (ar), ativa se dn/dT ≠ 0. Equivalente a tmhamtl-onoff.c: "
        "aproximação em que a superfície da amostra é livre; o fluido recebe calor e forma a TL, mas não "
        "esfria a amostra de volta. Para o acoplamento completo (Flizikowski et al., Opt. Express 28, "
        "7116, 2020) use TM+TL-Fluid, HAM (acoplamento completo): nos três metais do artigo a diferença é "
        "de 1e-4 a 1e-3 do sinal com o laser ligado e chega a ~1% logo após desligar (inox em ar).",
    "Thermal mirror of an opaque sample and thermal lens of the non-absorbing fluid in contact, with "
    "sample–fluid heat exchange (without the restriction k_f√D ≪ k√D_f). Lukasievicz et al., Appl. "
    "Spectrosc. 66, 1461 (2012), eqs. 20–22. k_s = D_s ρ_s c_s and k_f = D_f ρ_f c_f: fix ρ and c and "
    "fit the diffusivities. In air, the fluid barely changes the sample; in water, it does. It is the "
    "full model of Flizikowski et al., Opt. Express 28, 7116 (2020), eqs. 3–10; TM+TL-Fluid, HAM on/off "
    "(approximation) neglects the feedback of the fluid on the sample.":
        "Espelho térmico de amostra opaca e lente térmica do fluido não absorvente em contato, com troca "
        "de calor amostra–fluido (sem a restrição k_f√D ≪ k√D_f). Lukasievicz et al., Appl. Spectrosc. "
        "66, 1461 (2012), eqs. 20–22. k_s = D_s ρ_s c_s e k_f = D_f ρ_f c_f: fixe ρ e c e ajuste as "
        "difusividades. Para ar, o fluido quase não muda a amostra; em água, sim. É o modelo completo de "
        "Flizikowski et al., Opt. Express 28, 7116 (2020), eqs. 3–10; o TM+TL-Fluid, HAM on/off "
        "(aproximação) ignora o retorno do fluido sobre a amostra.",
    # ------------------------------------------------------------------ θ definitions
    "θ sample = ϑ·10⁻⁶,  ϑ = P_e A_e φ / (k_s λ_p)": "θ amostra = ϑ·10⁻⁶,  ϑ = P_e A_e φ / (k_s λ_p)",
    "θ fluid = 2 P_e A_e φ (∂n_f/∂T) / λ_p": "θ fluido = 2 P_e A_e φ (∂n_f/∂T) / λ_p",
    "The sample phase is ϑ·l·χ(α,l), with χ (dn/dT, α_T, ν, E, q∥, q⊥) in 10⁻⁶/K. It is not the θ in rad "
    "of the Shen model: the thickness l and the coefficient χ are left out of it.":
        "A fase da amostra é ϑ·l·χ(α,l), com χ (dn/dT, α_T, ν, E, q∥, q⊥) em 10⁻⁶/K. Não é o θ em rad do "
        "modelo de Shen: a espessura l e o coeficiente χ ficam fora dele.",
    "θ fluid / θ sample = 2·10⁶ k_s (∂n_f/∂T): with k_s known, it gives the fluid's ∂n/∂T.":
        "θ fluido / θ amostra = 2·10⁶ k_s (∂n_f/∂T): com k_s conhecido, dá o ∂n/∂T do fluido.",
    "θ of the classical thermal lens: includes the thickness L. Do not confuse it with the θ_TM of the "
    "thermal mirror (1/m) nor with the θ sample of TLLAM. ds/dT > 0: θ < 0 and the signal rises.":
        "θ da lente térmica clássica: inclui a espessura L. Não confundir com o θ_TM do espelho térmico "
        "(1/m) nem com o θ amostra do TLLAM. ds/dT > 0: θ < 0 e o sinal sobe.",
    "θ = −P_e β L φ (dn/dT) / (k λ_p),  β = absorption at t = 0 (before bleaching)":
        "θ = −P_e β L φ (dn/dT) / (k λ_p),  β = absorção em t = 0 (antes de branquear)",
    "Same θ as the Shen TL, but with the initial β: the long-time signal corresponds to ε̄·θ. dn/dT < 0 "
    "(water): θ > 0 and the signal drops. In a series, θ ∝ P_e. K_T [1/s] = K_m + K_pb P_e (K_m: mass "
    "diffusion; K_pb: photobleaching) and ε̄ = ε_products/ε_reactants — the K_T × P_e line appears in "
    "the series.":
        "Mesmo θ da TL de Shen, mas com o β inicial: o sinal em t longo corresponde a ε̄·θ. dn/dT < 0 "
        "(água): θ > 0 e o sinal cai. Na série, θ ∝ P_e. K_T [1/s] = K_m + K_pb P_e (K_m: difusão de "
        "massa; K_pb: fotobranqueamento) e ε̄ = ε_produtos/ε_reagentes — a reta de K_T × P_e aparece na "
        "série.",
    "Same θ as Shen. ds/dT < 0 (ZBLAN) with heating: θ_th > 0; laser cooling reverses the sign (φ < 0).":
        "Mesmo θ do Shen. ds/dT < 0 (ZBLAN) com aquecimento: θ_th > 0; resfriamento a laser inverte o "
        "sinal (φ < 0).",
    "Strength of the population lens. Δα_P: electronic polarizability difference excited–ground state; "
    "f_L = (n²+2)/3; E = h c/λ_ex; I0 = 2P/(π w_e²).":
        "Força da lente de população. Δα_P: diferença de polarizabilidade eletrônica excitado–fundamental; "
        "f_L = (n²+2)/3; E = h c/λ_ex; I0 = 2P/(π w_e²).",
    "Same θ_TM in the LAM, BAM and finite-sample models (it is the θ_TM of the papers, in m⁻¹). The HAM "
    "model uses θ_TM/A_e, dimensionless. Expanding sample (α_T > 0): θ_TM < 0.":
        "Mesmo θ_TM nos modelos LAM, BAM e amostra finita (é o θ_TM dos artigos, em m⁻¹). O modelo HAM usa "
        "θ_TM/A_e, adimensional. Amostra que dilata (α_T > 0): θ_TM < 0.",
    "Does not contain A_e (opaque sample). To compare with LAM, BAM and finite sample, which give θ_TM in "
    "1/m: θ_TM = θ_HAM·A_e. Expanding sample (α_T > 0): θ < 0. With incident P_e and reflectivity R, the "
    "θ_TM of Flizikowski et al. (2020) is (1−R)·θ_HAM.":
        "Não contém A_e (amostra opaca). Para comparar com LAM, BAM e amostra finita, que dão θ_TM em 1/m: "
        "θ_TM = θ_HAM·A_e. Amostra que dilata (α_T > 0): θ < 0. Com P_e incidente e refletividade R, o "
        "θ_TM de Flizikowski et al. (2020) é (1−R)·θ_HAM.",
    "Same θ as TM-HAM (opaque sample): does not contain A_e. Expanding sample: h_TM < 0. With incident "
    "P_e and reflectivity R, the θ_TM of Flizikowski et al. (2020) is (1−R)·h_TM.":
        "Mesmo θ do TM-HAM (amostra opaca): não contém A_e. Amostra que dilata: h_TM < 0. Com P_e incidente "
        "e refletividade R, o θ_TM de Flizikowski et al. (2020) é (1−R)·h_TM.",
    "Thermal lens of the fluid (two passes: before and after the reflection). With k_s known, h_f/P_e "
    "gives (∂n/∂T)_f: h_f = 2 P_e φ (∂n/∂T)_f/λ_p. P_e is the absorbed power; with incident P_e use "
    "P_e(1−R).":
        "Lente térmica do fluido (dois passes: antes e depois da reflexão). Com k_s conhecido, h_f/P_e dá "
        "(∂n/∂T)_f: h_f = 2 P_e φ (∂n/∂T)_f/λ_p. P_e é a potência absorvida; com P_e incidente use "
        "P_e(1−R).",
    # ------------------------------------------------------------------ fit log (optim.py)
    "equal bounds for {name}": "limites iguais para {name}",
    "none": "nenhum",
    "Two-level fit: sweep over {nl} ({n} points); {lin} fitted at each point":
        "Ajuste em dois níveis: varredura em {nl} ({n} pontos); {lin} ajustados em cada ponto",
    "Sweep: {n} local minima; refining {k} (χ² {chis})":
        "Varredura: {n} mínimos locais; refinando {k} (χ² {chis})",
    "Sweep: a single minimum (χ² {chi:.3e})": "Varredura: um único mínimo (χ² {chi:.3e})",
    "{n} free nonlinear parameters ({names}): using swarm (PSO) + LM":
        "{n} parâmetros não lineares livres ({names}): usando enxame (PSO) + LM",
    "Levenberg–Marquardt (normal accuracy, all parameters): χ² = {chi:.4e}":
        "Levenberg–Marquardt (precisão normal, todos os parâmetros): χ² = {chi:.4e}",
    "Run {run} of {n}: ": "Rodada {run} de {n}: ",
    "{tag}PSO with {particles} particles × {iters} iterations":
        "{tag}PSO com {particles} partículas × {iters} iterações",
    "{tag}PSO finished with χ² = {chi:.4e}": "{tag}PSO terminou com χ² = {chi:.4e}",
    "{tag}Levenberg–Marquardt finished with χ² = {chi:.4e}":
        "{tag}Levenberg–Marquardt terminou com χ² = {chi:.4e}",
    "All {n} runs reached the same minimum.": "As {n} rodadas chegaram ao mesmo mínimo.",
    "Warning: only {same} of {n} runs reached the best minimum (the others stopped at χ² = {chis}). "
    "There is more than one minimum; if the result looks odd, increase the number of runs or narrow "
    "the bounds.":
        "Atenção: só {same} de {n} rodadas chegaram ao melhor mínimo (as outras pararam em χ² = {chis}). "
        "Há mais de um mínimo; se o resultado parecer estranho, aumente o número de rodadas ou estreite "
        "os limites.",
    "Total: {n} model evaluations in {s:.1f} s": "Total: {n} avaliações do modelo em {s:.1f} s",
    # ------------------------------------------------------------------ series (batch.py)
    "{name}: failed ({err})": "{name}: falhou ({err})",
    "  (at bound: {names})": "  (no limite: {names})",
    "{k}/{n} {name}: {vals}; χ² = {chi:.3e} ({s:.1f} s){warn}":
        "{k}/{n} {name}: {vals}; χ² = {chi:.3e} ({s:.1f} s){warn}",
    "Series with {n} transients; free parameters: {free}":
        "Série com {n} transientes; parâmetros livres: {free}",
    "Series finished in {s:.1f} s": "Série concluída em {s:.1f} s",
    "Series stopped; the transients already fitted were kept.":
        "Série interrompida; os transientes já ajustados foram mantidos.",
    "ERROR: {err}": "ERRO: {err}",
    "transient {k}": "transiente {k}",
    # ------------------------------------------------------------------ server and data reading
    "Fit stopped by the user; keeping the best point found.":
        "Ajuste interrompido pelo usuário; mantendo o melhor ponto encontrado.",
    "example not found": "exemplo não encontrado",
    "unknown series": "série desconhecida",
    "unknown job": "job desconhecido",
    "not found": "não encontrado",
    "select at least one parameter to fit": "marque ao menos um parâmetro para ajustar",
    "no transient with at least 3 points": "nenhum transiente com pelo menos 3 pontos",
    "empty file": "arquivo vazio",
    "no numeric line found": "nenhuma linha numérica encontrada",
    "column {i}": "coluna {i}",
}
