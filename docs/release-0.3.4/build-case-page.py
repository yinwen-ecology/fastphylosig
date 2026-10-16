from pathlib import Path
import csv, json, re

root = Path(__file__).resolve().parents[2]
assets = root/'docs/assets/fast-hetero-case'
def rows(name):
    with (assets/name).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))
def number(value,digits=4):
    return 'NA' if value in {'NA',''} else f'{float(value):.{digits}g}'
def table(headers,data):
    return '\n'.join(['| '+' | '.join(headers)+' |','|'+'|'.join(['---']*len(headers))+'|']+
                     ['| '+' | '.join(str(x) for x in row)+' |' for row in data])

clades = {m:rows(m+'-clades.csv') for m in ['K','lambda','D','Delta']}
overall = {m:rows(m+'-overall.csv')[0] for m in clades}
hetero = {m:rows(m+'-heterogeneity.csv')[0] for m in clades}
labels = {
    'en': {
        'summary': ['Measure','Whole-tree estimate','Ordinary whole-tree signal P','Whole-tree n','Estimable / selected','Crown dispersion','Heterogeneity reference P'],
        'local': ['Crown','Original node','Continuous n','K','λ','D','Δ','Status'],
        'diagnostics': ['Crown','State counts','rank R-hat(Δ)','bulk ESS','tail ESS','Relative MCSE'],
        'descriptive': 'Not run (description)', 'single': 'D/Δ single-state',
        'estimable': 'Estimable', 'overall': 'Whole tree', 'states': 'Three states',
        'note': 'The K/λ reference P of 0.001 is the +1 rank result for 999 simulations with zero exceedances, and the resolution limit of this run. The ordinary whole-tree Δ P is 4/199 = 0.0201, with MCSE about 0.00995.',
        'page': 'fast-hetero-case.md',
    },
    'zh': {
        'summary': ['指标','整树估计','普通整树信号 P','整树 n','局部可用/选定','冠群离散','异质性参考 P'],
        'local': ['冠群','原始节点','连续性状 n','K','λ','D','Δ','状态'],
        'diagnostics': ['冠群','状态计数','rank R-hat(Δ)','bulk ESS','tail ESS','相对 MCSE'],
        'descriptive': '未计算（描述）', 'single': 'D/Δ 单状态',
        'estimable': '可估计', 'overall': '整体', 'states': '三状态',
        'note': 'K/λ 的参考 P = 0.001 是 999 次模拟、零次超越时的 +1 秩结果，也是这次模拟的分辨率下限。Δ 的普通整体 P 为 4/199 = 0.0201，MCSE 约为 0.00995。',
        'page': 'fast-hetero-case-zh.md',
    },
}
for lang, label in labels.items():
    summary = []
    for m in clades:
        g,h = overall[m],hetero[m]
        p = 'random: 0/1000; BM: '+number(g['P_Brownian']) if m=='D' else number(g['P'])
        summary.append([m,number(g['estimate']),p,g['n_species'],
          f'{h["n_clades"]} / 8',h['scale']+' = '+number(h['dispersion']),
          number(h['P_heterogeneity']) if m in ['K','lambda'] else label['descriptive']])
    summary_text = table(label['summary'],summary)+'\n\n'+label['note']
    local = []
    for i in range(8):
        local.append([f'C{i+1}',clades['K'][i]['node'],clades['K'][i]['n_species'],
          number(clades['K'][i]['estimate']),number(clades['lambda'][i]['estimate']),
          number(clades['D'][i]['estimate']),number(clades['Delta'][i]['estimate']),
          label['single'] if i==7 else label['estimable']])
    local_text = table(label['local'],local)
    diagnostics = []
    for z in clades['Delta']:
        if z['estimate']=='NA':
            continue
        diagnostics.append([z['crown'],z['state_counts'],number(z['rank_Rhat_Delta']),
          str(round(float(z['ESS_bulk_Delta']))),str(round(float(z['ESS_tail_Delta']))),
          f'{float(z["relative_MCSE_Delta"])*100:.2f}%'])
    g = overall['Delta']
    diagnostics.append([label['overall'],label['states'],number(g['rank_Rhat_Delta']),str(round(float(g['ESS_bulk_Delta']))),
                        str(round(float(g['ESS_tail_Delta']))),f'{float(g["relative_MCSE_Delta"])*100:.2f}%'])
    diagnostic_text = table(label['diagnostics'],diagnostics)
    page = root/'docs'/label['page']
    text = page.read_text(encoding='utf-8')
    for marker,value in [('SUMMARY',summary_text),('CLADES',local_text),('DIAGNOSTICS',diagnostic_text)]:
        pattern = f'<!-- GENERATED_{marker}_START -->.*?<!-- GENERATED_{marker}_END -->'
        text,count = re.subn(pattern,f'<!-- GENERATED_{marker}_START -->\n\n{value}\n\n<!-- GENERATED_{marker}_END -->',text,flags=re.S)
        assert count==1, (lang, marker)
    page.write_text(text,encoding='utf-8')
(assets/'case-summary.json').write_text(json.dumps({'version':'0.3.4','source':'native fast_hetero outputs',
    'methods':{m:{'overall':overall[m],'heterogeneity':hetero[m]} for m in clades}},ensure_ascii=False,indent=2),encoding='utf-8')
print('Generated English and Chinese tables from the same four saved result CSVs; no analysis rerun.')
