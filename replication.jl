include("ctRoot.jl")

function replicate_all(;saveall=false)

    fig1 = plot_announcement;(slides)
    saveall && savefig(fig1, "Graphs/current/announcements_paper.pdf", width = 900, height = 450)

    mt, dk = load("Output/mt.jld2", "mt", "dk")

    fig2 = Lplot(mt.ct; slides)
    saveall && savefig(fig2, "Graphs/current/first_L_paper.pdf", width = 900, height = 450)

    fig3 = hplot(mt.ct; slides)
    saveall && savefig(fig3, "Graphs/current/first_g_paper.pdf", width = 900, height = 350)

    fig4 = Eplot(mt.ct; slides)
    saveall && savefig(fig4, "Graphs/current/first_p_paper.pdf", width = 900, height = 400)

    fig5 = plansp(mt; slides)
    saveall && savefig(fig5, "Graphs/current/plans_paper.pdf", width = 900, height = 400)

    fig6 = plansp(mt; slides)
    saveall && savefig(fig6, "Graphs/current/planswide_paper.pdf", width=900, height=400)

    fig6 = Lωplot(mt; slides)
    saveall && savefig(fig6, "Graphs/current/contour_paper.pdf", width = 900, height = 500)

    fig7 = Cplot(mt; slides)
    saveall && savefig(fig7, "Graphs/current/Ccontour_paper.pdf", width = 900, height = 450)

    fig8 = avgplans(mt, 100, CIs = false; slides)
    saveall && savefig(fig8, "Graphs/current/mimics_paper.pdf", width = 900, height = 400)

    mt = load("Output/mt.jld2", "mt")

    fig9a, fig9b = strategy_μ(mt; slides, save_stats = true, folder = "Graphs/current/");
    saveall && savefig(fig9a, "Graphs/current/marg_achi_paper.pdf", width = 600, height = 450)
    saveall && savefig(fig9b, "Graphs/current/marg_omegachi_paper.pdf", width = 600, height = 450)

    # Comparative statics would go here

    fig12 = comp_plot_planner(mt; slides, k = 5, Average = true, Kambe = true)
    saveall && savefig(fig12, "Graphs/current/ramsey_paper.pdf", width = 900, height = 350)

    fig_app1 = Lplot_fixed_ω(mt; slides)
    saveall && savefig(fig_app1, "Graphs/current/contour_app_paper.pdf", width = 900, height = 400)


    # Gradual Feedback Rules

    mtp = load("Output/mtpsi_s20_20230720.jld2", "mt_psi");

    fig13a = twolines(mtp, jp = 2; slides)
    saveall && savefig(fig13a, "Graphs/current/psithing_L1.pdf", width = 500, height = 350)

    fig13b = twolines(mtp, show = "C", jp = 2; slides)
    saveall && savefig(fig13b, "Graphs/current/psithing_C1.pdf", width = 500, height = 350)

    fig14 = implied_plan(mtp; slides)
    saveall && savefig(fig14, "Graphs/current/psithing_plans.pdf", width = 1000, height = 450)
end

function thumbnail(;saveall = false)

    mt = load("Output/mt.jld2", "mt")

    fig = Lplot(mt.ct, slides = true, plot_bgcolor="#fafafa", paper_bgcolor="#fafafa")
    saveall && savefig(fig, "cddp_fig.png", width = 850, height = 500)
end


function make_slides(; folder = "Graphs/current/", saveall = false)
    slides = true
    qual = "_slides"

    fig1 = plot_announcements(; slides)
    saveall && savefig(fig1, folder * "/announcements"*qual*".pdf", width = 600, height = 450)

    fig1b = plot_announcements(; slides, cond = true)
    saveall && savefig(fig1b, folder * "/announcements_cond"*qual*".pdf", width = 600, height = 450)

    fig1c = plot_announcements(; slides, cond_t = true)
    saveall && savefig(fig1c, folder * "/announcements_cond_t"*qual*".pdf", width = 600, height = 450)

    mt = load("Output/mt.jld2", "mt")

    fig2 = Lplot(mt.ct; slides)
    saveall && savefig(fig2, folder * "/first_L"*qual*".pdf", width = 800, height = 550)

    fig3 = hplot(mt.ct; slides)
    saveall && savefig(fig3, folder * "/first_g"*qual*".pdf", width = 900, height = 550)

    fig4 = Eplot(mt.ct; slides)
    saveall && savefig(fig4, folder * "/first_p"*qual*".pdf", width = 900, height = 550)

    fig7 = Cplot(mt; slides)
    saveall && savefig(fig7, folder * "Ccontour"*qual*".pdf", width = 900, height = 550)
    
    fig5 = plansp(mt; slides)
    saveall && savefig(fig5, folder * "plans"*qual*".pdf", width = 900, height = 600)

    fig6 = Lωplot(mt; slides)
    saveall && savefig(fig6, folder * "contour"*qual*".pdf", width = 900, height = 550)

    fig9a, fig9b = strategy_μ(mt; slides, save_stats = true, folder);
    saveall && savefig(fig9a, folder * "marg_achi"*qual*".pdf", width = 600, height = 450)
    saveall && savefig(fig9b, folder * "marg_omegachi"*qual*".pdf", width = 600, height = 450)


    figRamsey = comp_plot_planner(mt; slides, k=5, Average=false, Kambe=false, title = "", share = true)
    saveall && savefig(figRamsey, folder * "Ramseyplan" * qual * ".pdf", width=900, height=550)

    figRamseyKambe = comp_plot_planner(mt; slides, Kambe = true, Average = false)
    saveall && savefig(figRamseyKambe, folder * "Ramsey_optimal" * qual * ".pdf", width = 900, height = 550)

end


function replicate_JET2(; saveall = false)
    folder = "Graphs/JET/"
    qual = "_paper"
    slides = false
    mt, dk = load("Output/JET/mt.jld2", "mt", "dk")

    mtp = load("Output/JET/mtpsi", "mt_psi")

    figRamsey = comp_plot_planner(mt; slides, k=5, Average=false, Kambe=false, title="", share=true)
    saveall && savefig(figRamsey, folder * "Ramseyplan" * qual * ".pdf", width=850, height=350)

    fig2 = plot_announcements(;slides)
    saveall && savefig(fig2, folder * "announcements" * qual * ".pdf", width=900, height=450)


    fig3 = hplot(mt.ct; slides)
    saveall && savefig(fig3, folder * "first_g" * qual * ".pdf", width=900, height=350)

    fig4 = Lplot(mt.ct; slides)
    saveall && savefig(fig4, folder * "first_L" * qual * ".pdf", width=900, height=450)

    fig5 = planspwide(mt, dk; slides, share = true)
    saveall && savefig(fig5, folder * "planswide" * qual * ".pdf", width=900, height=400)

    fig6 = Lωplot(mt; slides)
    saveall && savefig(fig6, folder * "contour" * qual * ".pdf", width=900, height=500)

    fig7 = Cplot(mt; slides)
    saveall && savefig(fig7, folder * "Ccontour" * qual * ".pdf", width=900, height=450)


    # Comparative statics
    fig8_old = plot_cs(:σ; slides)
    saveall && savefig(fig8_old, folder * "compstats_sigma_old" * qual * ".pdf", width=900, height=400)

    fig9b_old = plot_cs(:κ; slides)
    saveall && savefig(fig9b_old, folder * "compstats_kappa_old" * qual * ".pdf", width=500, height=350)

    # fig8 = plot_compstats(:σ)
    # saveall && savefig(fig8, folder * "compstats_sigma" * qual * ".pdf", width = 900, height = 400)

    fig13_app = Eplot(mt.ct; slides)
    saveall && savefig(fig13_app, folder * "first_p" * qual * ".pdf", width=900, height=400)

    fig14_app = planspwide(mt; slides, share=true)
    saveall && savefig(fig14_app, folder * "planswide_r" * qual * ".pdf", width=900, height=400)


    fig13a = twolines(mtp, jp=2; slides)
    saveall && savefig(fig13a, folder * "psithing_L1" * qual* ".pdf", width=500, height=350)

    fig13b = twolines(mtp, show="C", jp=2; slides)
    saveall && savefig(fig13b, folder * "psithing_C1" * qual* ".pdf", width=500, height=350)

    fig14 = implied_plan_wide(mtp; slides)
    saveall && savefig(fig14, folder * "psithing_plans" * qual* ".pdf", width=900, height=350)

    fig16_old = plot_cs(:σ; showL = true, slides)
    saveall && savefig(fig16_old, folder * "compstats_L_sigma_old" * qual * ".pdf", width=900, height=400)

end
