"""Render the bilingual fast_hetero workflow as self-contained, editable SVGs."""
from pathlib import Path
from html import escape

root = Path(__file__).resolve().parents[2] / 'docs'
copy = {
    'en': {
        'title': 'fast_hetero() · Technical workflow',
        'subtitle': 'From a shared tree to local signal, whole-tree context and biological interpretation',
        'main': 'Local descriptive path', 'optional': 'Optional analysis',
        'input': ['01 · Input and compile the tree', 'Rooted tree with branch lengths + species-named traits; choose K, λ, D or Δ.', 'Compile topology and branch geometry once for the crown scan.'],
        'select': ['02 · Select crown clades', 'Prespecified nodes, or maximal disjoint crowns', 'within min_tips / max_tips bounds.', 'Retain node IDs, coverage and overlap.'],
        'match': ['03 · Align observed traits', 'Match species; remove missing values by trait.', 'Reuse workspaces for shared missing-data masks.', 'Check observed counts and categorical states.'],
        'kernel': '04 · Scan crowns with shared C++ kernels',
        'cards': [['K · continuous', 'fast_k kernels', 'Contrasts / permutations'], ['λ · continuous', 'fast_lambda kernels', 'Likelihood / optimization'], ['D · binary', 'fast_d kernels', 'Contrasts / calibrations'], ['Δ · categorical', 'fast_ace → node P', 'fast_delta → MCMC']],
        'overall': ['Optional A · Whole-tree analysis', 'overall = TRUE: estimate on the full tree.', 'overall_test = TRUE: add its ordinary signal test.', 'Reuse the corresponding method kernels.'],
        'reference': ['Optional B · Dispersion reference', 'test = TRUE; use disjoint crowns.', 'K: BM · λ: approximate common-local-λ.', 'D: count-conditioned ER (experimental).'],
        'output': ['05 · Summarize and plot', 'clades + heterogeneity + selection; overall when requested.', 'Local/overall contrasts, state composition, fit diagnostics and failures remain visible.', 'Native plot(): crown colours, separate test annotations; plotting preserves the RNG state.'],
        'footer': 'Recommended Δ workflow: local and whole-tree description, with fit/MCMC diagnostics.',
    },
    'zh': {
        'title': 'fast_hetero() · 技术路线',
        'subtitle': '从共享系统发育树到局部信号、总体背景与生物学解释',
        'main': '局部描述主流程', 'optional': '按需启用的分析',
        'input': ['01 · 输入与树编译', '带枝长的有根树 + 按物种命名的性状；指定 K、λ、D 或 Δ。', '为冠群扫描一次性编译拓扑与枝长结构。'],
        'select': ['02 · 选择单系冠群', '指定节点，或在 min_tips / max_tips 范围内', '自动选择最大不重叠冠群。', '保留原始节点、覆盖度与重叠信息。'],
        'match': ['03 · 对齐并筛选实际观测', '匹配物种，按性状处理缺失值。', '对相同缺失模式复用工作空间。', '检查实际观测数与分类状态组成。'],
        'kernel': '04 · 共享 C++ 内核扫描冠群',
        'cards': [['K · 连续性状', 'fast_k 内核', '对比 / 置换'], ['λ · 连续性状', 'fast_lambda 内核', '似然 / 优化'], ['D · 二元性状', 'fast_d 内核', '对比 / 参考模拟'], ['Δ · 分类性状', 'fast_ace → 祖先概率', 'fast_delta → MCMC']],
        'overall': ['可选 A · 整树分析', 'overall = TRUE：计算整树信号估计。', 'overall_test = TRUE：增加普通整树信号检验。', '复用对应指标的计算内核。'],
        'reference': ['可选 B · 冠群离散参考', 'test = TRUE；使用不重叠冠群。', 'K：BM；λ：共同局部 λ 的近似参考。', 'D：计数条件 ER（实验性）。'],
        'output': ['05 · 汇总与绘图', 'clades + heterogeneity + selection；按需返回 overall。', '保留局部与总体差异、状态组成、拟合诊断和失败记录。', '原生 plot()：冠群着色、检验分别标注；绘图保持随机数状态。'],
        'footer': 'Δ 推荐流程：局部与总体描述，同时检查拟合和 MCMC 诊断。',
    },
}

def render(lang, c):
    parts = [f'''<svg xmlns="http://www.w3.org/2000/svg" width="1120" height="1080" viewBox="0 0 1120 1080" role="img" aria-labelledby="title desc" xml:lang="{'zh-CN' if lang == 'zh' else 'en'}">
<title id="title">{escape(c['title'])}</title><desc id="desc">{escape(c['subtitle'])}. {escape(c['footer'])}</desc>
<defs><marker id="arrow" markerWidth="8" markerHeight="8" refX="7" refY="4" orient="auto"><path d="M0 0 L8 4 L0 8Z" fill="#3e7888"/></marker></defs>
<style>text{{font-family:system-ui,-apple-system,'Segoe UI','Microsoft YaHei',sans-serif;fill:#17364d}} .title{{font-size:32px;font-weight:750}} .heading{{font-size:23px;font-weight:700}} .body{{font-size:18px}} .muted{{font-size:17px;fill:#52687a}} .card-title{{font-size:20px;font-weight:700}} .card-text{{font-size:17px}} .flow{{fill:none;stroke:#3e7888;stroke-width:2.2;marker-end:url(#arrow)}} .optional{{stroke-dasharray:7 5}}</style>
<rect width="1120" height="1080" rx="16" fill="#fff"/>
''']
    def text(x, y, value, cls='body'):
        parts.append(f'<text x="{x}" y="{y}" class="{cls}">{escape(value)}</text>')
    def box(x, y, w, h, fill='#f0f7f9', dashed=False):
        parts.append(f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="12" fill="{fill}" stroke="#a9c7d0" stroke-width="1.3"'+(' stroke-dasharray="7 5"' if dashed else '')+'/>')
    def flow(d, optional=False):
        parts.append(f'<path d="{d}" class="flow'+(' optional' if optional else '')+'"/>')
    def lines(x,y,values):
        text(x,y,values[0],'heading')
        for i,value in enumerate(values[1:]):text(x,y+32+27*i,value)
    text(40,47,c['title'],'title');text(40,79,c['subtitle'],'muted')
    flow('M42 105 H83');text(95,110,c['main'],'muted')
    flow('M365 105 H406',True);text(418,110,c['optional'],'muted')
    box(40,135,1040,100);lines(65,167,c['input'])
    flow('M290 235 V265')
    box(40,265,500,135);lines(65,298,c['select'])
    flow('M540 332 H580')
    box(580,265,500,135);lines(605,298,c['match'])
    flow('M830 400 V417 H560 V440')
    box(40,440,1040,155,fill='#edf6f1');text(65,475,c['kernel'],'heading')
    for i, card in enumerate(c['cards']):
        x=55+258*i
        box(x,491,235,88,fill='#fff')
        text(x+14,516,card[0],'card-title');text(x+14,542,card[1],'card-text');text(x+14,565,card[2],'card-text')
    # Keep a direct descriptive path between the two optional panels.
    flow('M560 595 V855')
    flow('M290 595 V640',True);flow('M830 595 V640',True)
    box(40,640,500,160,fill='#f8fafc',dashed=True);lines(65,675,c['overall'])
    box(580,640,500,160,fill='#f8fafc',dashed=True);lines(605,675,c['reference'])
    flow('M290 800 V855',True);flow('M830 800 V855',True)
    box(40,855,1040,145);lines(65,890,c['output'])
    text(40,1040,c['footer'],'muted')
    parts.append('</svg>\n')
    return ''.join(parts)

for lang,c in copy.items():
    path=root/('fast-hetero-roadmap'+('-zh' if lang=='zh' else '')+'.svg')
    path.write_text(render(lang,c),encoding='utf-8')
print('Rendered English and Chinese fast_hetero SVG workflows.')
