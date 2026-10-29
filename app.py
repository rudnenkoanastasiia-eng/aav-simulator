import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

st.set_page_config(
    page_title="Навчальний симулятор rAAV", page_icon="🧬", layout="wide"
)

st.title("🧬 Професійний симулятор біопроцесу rAAV (HEK293)")
st.caption(
    "Модель upstream- та downstream-виробництва вектора rAAV методом потрійної"
    " трансфекції."
)

# --- ТАБИ НАВІГАЦІЇ ---
tab1, tab2, tab3, tab4 = st.tabs([
    "🎮 Симулятор та розрахунки",
    "🖼️ Обладнання та GMP-лінія",
    "📚 Картки-шпори (Покращення якості)",
    "⚖️ Масштабування (200L ➔ 2000L)",
])

# --- SIDEBAR: ВХІДНІ ПАРАМЕТРИ ---
st.sidebar.header("⚙️ 1. Параметри клітин та ДНК")

cell_density = st.sidebar.slider(
    "Щільність клітин при трансфекції (Xₜ), 10⁶ cells/mL",
    0.5,
    8.0,
    2.0,
    0.1,
    help="Суспензійний оптимум HEK293: 1.5–2.5 x 10⁶ cells/mL",
)

dna_per_cell = st.sidebar.slider(
    "Концентрація ДНК, мкг / 10⁶ клітин",
    0.2,
    3.0,
    1.2,
    0.1,
    help="Оптимум: 1.0–1.5 мкг на мільйон клітин",
)

pei_ratio = st.sidebar.slider(
    "Співвідношення ДНК : PEI (масове w/w)",
    1.0,
    5.0,
    3.0,
    0.1,
    help="Оптимум: від 1:2 до 1:4 (w/w)",
)

complex_time = st.sidebar.slider(
    "Час формування комплексу DNA:PEI, хв",
    5,
    45,
    18,
    1,
    help="Оптимум: 15–20 хвилин. Понад 25 хв викликає агломерацію",
)

st.sidebar.markdown("**Співвідношення плазмід (Сума = 1.0)**")
ratio_transgene = st.sidebar.slider(
    "Касета трансгена (GoI)", 0.1, 0.7, 0.33, 0.01
)
ratio_repcap = st.sidebar.slider("Rep/Cap плазміда", 0.1, 0.7, 0.33, 0.01)
ratio_helper = st.sidebar.slider("Helper плазміда", 0.1, 0.7, 0.34, 0.01)

st.sidebar.header("🌡️ 2. Умови культивування")

temp_mode = st.sidebar.selectbox(
    "Температурний режим після трансфекції",
    options=[
        "Постійна температура (37.0 °C)",
        "Гіпотермічний шифт (33.0 °C після трансфекції)",
    ],
    help="Зниження до 33°C покращує згортання білків та збірку капсидів",
)

harvest_hpt = st.sidebar.slider(
    "Час збору (Harvest time), hpt",
    24,
    120,
    60,
    1,
    help="Оптимальний збір: 48–72 години після трансфекції",
)

serotype = st.sidebar.selectbox(
    "Серотип rAAV", options=["AAV2", "AAV5", "AAV8", "AAV9"]
)

lysis_mode = st.sidebar.selectbox(
    "Стратегія лізису клітин",
    options=[
        "Без лізису (збір супернатанту)",
        "М'який / Частковий лізис",
        "Повний лізис",
    ],
)

st.sidebar.header("🏭 3. Масштаб та Очищення")

scale_l = st.sidebar.select_slider(
    "Масштаб біореактора, L", options=[2, 50, 200, 500, 1000, 2000], value=200
)

purification_mode = st.sidebar.selectbox(
    "Метод Downstream-очищення",
    options=[
        "Без очищення (Crude Harvest / Сирий лізат) [Вихід 100%, Низька чистота]",
        "Ультрацентрифугування (Iodixanol/CsCl) [Вихід ~40%, Чистота ~85%]",
        "Хроматографія (Affinity + IEX) [Вихід ~65%, Чистота >99%]",
    ],
)

# --- РОЗРАХУНКОВИЙ БЛОК ---
base_titer = 3.5e10

cell_effect = max(0.1, 1.0 - abs(cell_density - 2.0) / 6.0)
dna_effect = max(
    0.1, 1.0 + 0.2 * (dna_per_cell - 1.2) - 0.25 * max(0.0, dna_per_cell - 1.5)
)
pei_effect = max(0.1, 1.0 - abs(pei_ratio - 3.0) * 0.2)

if 15 <= complex_time <= 20:
  complex_effect = 1.0
elif complex_time < 15:
  complex_effect = 0.7 + (complex_time / 15.0) * 0.3
else:
  complex_effect = max(0.2, 1.0 - (complex_time - 20) * 0.035)

plasmid_sum = ratio_transgene + ratio_repcap + ratio_helper
imbalance = abs(plasmid_sum - 1.0) + abs(ratio_transgene - 0.33) * 0.5
balance_effect = max(0.1, 1.0 - 1.5 * imbalance)

if harvest_hpt < 48:
  time_effect = 0.6 + (harvest_hpt - 24) / 24 * 0.4
elif harvest_hpt <= 72:
  time_effect = 1.0
else:
  time_effect = 1.0 - (harvest_hpt - 72) / 48 * 0.25

intracellular_fraction = {
    "AAV2": 0.85,
    "AAV5": 0.60,
    "AAV8": 0.75,
    "AAV9": 0.70,
}[serotype]

lysis_release = {
    "Без лізису (збір супернатанту)": 0.0,
    "М'який / Частковий лізис": 0.60,
    "Повний лізис": 0.95,
}[lysis_mode]

available_fraction = (1.0 - intracellular_fraction) + (
    intracellular_fraction * lysis_release
)
scale_penalty = {2: 1.00, 50: 0.98, 200: 0.95, 500: 0.90, 1000: 0.85, 2000: 0.78}[
    scale_l
]

titer = (
    base_titer
    * cell_effect
    * dna_effect
    * pei_effect
    * complex_effect
    * balance_effect
    * time_effect
    * available_fraction
    * scale_penalty
)

temp_bonus = 4.5 if temp_mode.startswith("Гіпотермічний") else 0.0
base_full = 0.28
full_fraction = (
    base_full
    - 0.07 * max(0.0, dna_per_cell - 1.2)
    - 0.15 * imbalance
    - 0.05 * max(0.0, harvest_hpt - 72) / 24.0
    - 0.03 * max(0.0, cell_density - 4.0)
)
full_pct = max(5.0, min(50.0, (full_fraction * 100) + temp_bonus))

impurity_load = (
    1.0
    + 2.0 * lysis_release
    + 0.5 * max(0.0, harvest_hpt - 72) / 24.0
    + 0.3 * max(0.0, cell_density - 4.0)
)

if purification_mode.startswith("Без очищення"):
  downstream_recovery = 1.00
  purification_penalty = 35.0
elif purification_mode.startswith("Ультрацентрифугування"):
  downstream_recovery = 0.40
  purification_penalty = 12.0
else:
  downstream_recovery = 0.65
  purification_penalty = 0.0

quality_score = (
    0.50 * (full_pct / 35.0 * 100)
    + 0.30 * max(0.0, (1.0 - (impurity_load - 1.0) / 3.0) * 100)
    + 0.20 * min(100.0, (titer / 3.5e10) * 100)
)
quality_score = max(0.0, min(100.0, quality_score - purification_penalty))

harvest_volume_ml = scale_l * 1000.0 * 0.9
crude_total_vg = titer * harvest_volume_ml
purified_total_vg = crude_total_vg * downstream_recovery

clean_purif_name = purification_mode.split("[")[0].strip()

# --- АНАЛІТИЧНИЙ МОДУЛЬ ДИНАМІЧНИХ РЕКОМЕНДАЦІЙ ---
recommendations = []

if not temp_mode.startswith("Гіпотермічний"):
  recommendations.append(
      "❄️ **Увімкніть Гіпотермічний шифт (33.0 °C):** Це підніме частку повних"
      " капсидів (% full) на **+4.5%** за рахунок покращення фолдингу білків та"
      " зменшення метаболічного стресу."
  )

if abs(plasmid_sum - 1.0) > 0.02 or abs(ratio_transgene - 0.33) > 0.05:
  recommendations.append(
      "⚖️ **Збалансуйте плазміди (1:1:1):** Установіть частку трансгена = 0.33,"
      " Rep/Cap = 0.33, Helper = 0.34. Поточний дисбаланс знижує частку повних"
      " часток та титр."
  )

if complex_time < 15 or complex_time > 20:
  recommendations.append(
      f"⏱️ **Скоригуйте час комплексу DNA:PEI:** Оптимум — **15–20 хвилин**"
      f" (зараз {complex_time} хв). { 'Час занадто малий для формування комплексу.' if complex_time < 15 else 'Понад 22 хв викликає агломерацію та втрату трансфекційності.' }"
  )

if dna_per_cell > 1.3:
  recommendations.append(
      f"🧬 **Зменшіть дозу ДНК:** Оптимум становить **1.0–1.2 мкг / 10⁶ клітин**"
      f" (зараз {dna_per_cell} мкг). Передозування ДНК виснажує клітини та"
      " знижує % full."
  )

if harvest_hpt > 72:
  recommendations.append(
      f"⏰ **Оптимізуйте час збору (Harvest):** Скоротіть час до **48–72 hpt**"
      f" (зараз {harvest_hpt} hpt), щоб уникнути масового відмирання клітин та"
      " накопичення домішок (HCDNA/HCP)."
  )

if purification_mode.startswith("Без очищення"):
  recommendations.append(
      "🧪 **Переключіть на Хроматографію (Affinity + IEX):** Сирий лізат має"
      " критично низький рівень чистоти і не може використовуватися in vivo"
      " через токсичність."
  )

if scale_l >= 1000:
  recommendations.append(
      "🏭 **Для масштабу 1000–2000 L:** Використовуйте автоматичне inline-змішування"
      " PEI:ДНК та оптимізуйте kLa/PV для захисту клітин від зсувного"
      " напруження (shear stress)."
  )

# --- ТАБ 1: СИМУЛЯТОР ТА РОЗРАХУНКИ ---
with tab1:
  st.subheader("📊 Головні показники лоту rAAV")
  c1, c2, c3, c4 = st.columns(4)
  c1.metric("Титр (Crude Harvest)", f"{titer:.2e} vg/mL")
  c2.metric("Вихід продукту (Purified/Total Yield)", f"{purified_total_vg:.2e} vg")
  c3.metric("Частка повних капсидів", f"{full_pct:.1f} %")
  c4.metric("Індекс якості лоту", f"{quality_score:.0f} / 100")

  st.markdown("---")
  st.subheader("⚠️ Попередження та критичні компроміси")

  if purification_mode.startswith("Без очищення"):
    st.error(
        "🚨 **Увага (Без очищення):** Сирий лізат містить величезну кількість"
        " клітинного сміття, залишків ДНК (HCDNA) та токсичних білків (HCP)."
        " Препарат непридатний для використання in vivo без хроматографії!"
    )

  if abs(plasmid_sum - 1.0) > 0.02:
    st.error(f"❌ Сума часток плазмід має дорівнювати 1.00 (Зараз: {plasmid_sum:.2f})")

  if complex_time > 22:
    st.error(
        f"⏱️ **Час формування комплексу DNA:PEI занадто великий ({complex_time}"
        " хв)!** Комплекси починають агломерувати, знижуючи ефективність"
        " трансфекції."
    )

  if pei_ratio < 1.8 or pei_ratio > 4.2:
    st.warning(
        f"⚠️ **Неоптимальне співвідношення DNA:PEI ({pei_ratio}:1):** Нестача PEI"
        " погіршує доставку ДНК, а надлишок є цитотоксичним."
    )

  if titer > 3.0e10 and full_pct < 20.0:
    st.warning(
        "⚡ **Високий титр, але низька якість!** Багато часток є порожніми (без"
        " генома). Геномний титр ≠ інфекційний титр."
    )

  if harvest_hpt > 84:
    st.warning(
        "⏱️ **Пізній збір (>84 hpt):** пришвидшує деградацію вірусу та масове"
        " накопичення ДНК/білків хазяйки (HCDNA/HCP)."
    )

  if temp_mode.startswith("Гіпотермічний"):
    st.success(
        "❄️ **Гіпотермічний шифт (33°C) активовано:** Покращено згортання білків,"
        " частку повних капсидів (+4.5%) та знижено метаболічний стрес клітин."
    )

  st.subheader("📈 Динаміка біопроцесу (0 ➔ Harvest)")
  t_steps = np.linspace(0, harvest_hpt, 100)
  xv_t = cell_density * np.exp(0.015 * t_steps) * (1 - t_steps / 150)
  vg_t = (
      (titer / (1 + np.exp(-0.1 * (t_steps - 36))))
      if harvest_hpt > 24
      else np.zeros(100)
  )

  fig = go.Figure()
  fig.add_trace(
      go.Scatter(
          x=t_steps,
          y=xv_t,
          name="Щільність клітин (10⁶/mL)",
          line=dict(color="blue", width=3),
      )
  )
  fig.add_trace(
      go.Scatter(
          x=t_steps,
          y=vg_t,
          name="Титр AAV (vg/mL)",
          yaxis="y2",
          line=dict(color="green", width=3, dash="dash"),
      )
  )

  fig.update_layout(
      title="Кінетика росту клітин та накопичення rAAV у біореакторі",
      xaxis=dict(title="Час після трансфекції (hpt), годин"),
      yaxis=dict(title="Живі клітини (10⁶ cells/mL)", side="left"),
      yaxis2=dict(
          title="Титр rAAV (vg/mL)", overlaying="y", side="right", type="log"
      ),
      legend=dict(orientation="h", y=1.15),
      height=400,
  )
  st.plotly_chart(fig, use_container_width=True)

# --- ТАБ 2: ВІЗУАЛІЗАЦІЯ БІОРЕАКТОРА ТА ОБЛАДНАННЯ GMP-ЛІНІЇ ---
with tab2:
  st.subheader("🏭 Графічне обладнання технологічного маршруту (Equipment PFD)")
  st.caption(
      "Схема ключових апаратів GMP-лінії з прямим відображенням заданих"
      " параметрів:"
  )

  eq1, eq2, eq3, eq4, eq5 = st.columns(5)

  target_temp = "33.0 °C" if temp_mode.startswith("Гіпотермічний") else "37.0 °C"

  # 1. SEED TRAIN (Шейкерна колба + Wave-реактор)
  with eq1:
    st.markdown("##### 1. Посівний етап\n**(Seed Train / Wave)**")
    svg_seed = """
        <svg width="240" height="290" viewBox="0 0 240 290" xmlns="http://www.w3.org/2000/svg">
          <rect x="5" y="5" width="230" height="280" rx="12" fill="#f8f9fa" stroke="#1a73e8" stroke-width="2.5"/>
          <text x="120" y="32" font-family="Arial" font-size="13" font-weight="bold" text-anchor="middle" fill="#1a73e8">Колби & Wave 10L</text>
          <text x="120" y="52" font-family="Arial" font-size="11" text-anchor="middle" fill="#202124">WCB HEK293 (Інокулят)</text>
          <text x="120" y="70" font-family="Arial" font-size="11" text-anchor="middle" fill="#3c4043">X₀ = 0.3×10⁶/mL</text>

          <!-- Erlenmeyer Flask -->
          <polygon points="60,130 35,190 95,190 70,130" fill="#e8f0fe" stroke="#1a73e8" stroke-width="3"/>
          <rect x="61" y="112" width="18" height="18" fill="#e8f0fe" stroke="#1a73e8" stroke-width="2"/>
          <path d="M38,183 Q65,173 92,183 L95,190 L35,190 Z" fill="#4285f4" opacity="0.7"/>

          <!-- Wave Bioreactor Bag -->
          <rect x="120" y="140" width="85" height="50" rx="10" fill="#e8f0fe" stroke="#1a73e8" stroke-width="3"/>
          <rect x="125" y="152" width="75" height="32" rx="6" fill="#4285f4" opacity="0.6"/>
          <line x1="110" y1="198" x2="215" y2="198" stroke="#5f6368" stroke-width="5"/>

          <text x="120" y="250" font-family="Arial" font-size="12" font-weight="bold" text-anchor="middle" fill="#137333">Життєздатність >98%</text>
        </svg>
        """
    st.components.v1.html(svg_seed, height=300)

  # 2. MAIN BIOREACTOR (Single-Use Stirred Tank)
  with eq2:
    st.markdown(f"##### 2. Біореактор SUB\n**({scale_l} L Production)**")
    shear_warning = scale_l >= 1000
    reactor_fill = "#e3f2fd" if not shear_warning else "#fff3e0"
    svg_sub = f"""
        <svg width="240" height="290" viewBox="0 0 240 290" xmlns="http://www.w3.org/2000/svg">
          <rect x="5" y="5" width="230" height="280" rx="12" fill="#f8f9fa" stroke="#137333" stroke-width="2.5"/>
          <text x="120" y="32" font-family="Arial" font-size="13" font-weight="bold" text-anchor="middle" fill="#137333">SUB {scale_l} L</text>

          <!-- Vessel -->
          <rect x="50" y="80" width="140" height="170" rx="18" fill="{reactor_fill}" stroke="#137333" stroke-width="3"/>
          <!-- Agitator Shaft & Motor -->
          <rect x="100" y="48" width="40" height="18" fill="#202124" rx="4"/>
          <line x1="120" y1="66" x2="120" y2="225" stroke="#3c4043" stroke-width="4"/>
          <rect x="85" y="190" width="70" height="12" fill="#5f6368" rx="3"/>
          <rect x="92" y="212" width="56" height="10" fill="#5f6368" rx="3"/>

          <!-- Readable Badge Overlay -->
          <rect x="65" y="100" width="110" height="68" rx="6" fill="#ffffff" opacity="0.9"/>
          <text x="120" y="118" font-family="Arial" font-size="11" text-anchor="middle" font-weight="bold" fill="#0d47a1">Xₜ: {cell_density}M/mL</text>
          <text x="120" y="138" font-family="Arial" font-size="11" text-anchor="middle" font-weight="bold" fill="#b71c1c">T: {target_temp}</text>
          <text x="120" y="156" font-family="Arial" font-size="10" text-anchor="middle" fill="#3c4043">PEI: {complex_time} хв</text>
        </svg>
        """
    st.components.v1.html(svg_sub, height=300)

  # 3. CLARIFICATION & LYSIS (Глибинний фільтр)
  with eq3:
    st.markdown("##### 3. Освітлення\n**(Depth Filter Capsule)**")
    svg_filter = f"""
        <svg width="240" height="290" viewBox="0 0 240 290" xmlns="http://www.w3.org/2000/svg">
          <rect x="5" y="5" width="230" height="280" rx="12" fill="#f8f9fa" stroke="#b06000" stroke-width="2.5"/>
          <text x="120" y="32" font-family="Arial" font-size="13" font-weight="bold" text-anchor="middle" fill="#b06000">Depth Filter 0.22µm</text>

          <!-- Filter Housing -->
          <rect x="65" y="75" width="110" height="160" rx="10" fill="#fef7e0" stroke="#b06000" stroke-width="3"/>
          <line x1="70" y1="115" x2="170" y2="115" stroke="#b06000" stroke-width="3" stroke-dasharray="5,3"/>
          <line x1="70" y1="150" x2="170" y2="150" stroke="#b06000" stroke-width="3" stroke-dasharray="5,3"/>
          <line x1="70" y1="185" x2="170" y2="185" stroke="#b06000" stroke-width="3" stroke-dasharray="5,3"/>

          <!-- Pipes -->
          <line x1="120" y1="50" x2="120" y2="75" stroke="#b06000" stroke-width="5"/>
          <line x1="120" y1="235" x2="120" y2="260" stroke="#b06000" stroke-width="5"/>

          <text x="120" y="136" font-family="Arial" font-size="11" font-weight="bold" text-anchor="middle" fill="#202124">{lysis_mode.split('(')[0]}</text>
          <text x="120" y="210" font-family="Arial" font-size="12" font-weight="bold" text-anchor="middle" fill="#b71c1c">HCDNA: {impurity_load:.1f}x</text>
        </svg>
        """
    st.components.v1.html(svg_filter, height=300)

  # 4. CHROMATOGRAPHY SKID (Афінна та Аніонообмінна колонка)
  with eq4:
    st.markdown("##### 4. Хроматографія\n**(Affinity & AEX Column)**")
    svg_chrom = f"""
        <svg width="240" height="290" viewBox="0 0 240 290" xmlns="http://www.w3.org/2000/svg">
          <rect x="5" y="5" width="230" height="280" rx="12" fill="#f8f9fa" stroke="#c5221f" stroke-width="2.5"/>
          <text x="120" y="32" font-family="Arial" font-size="13" font-weight="bold" text-anchor="middle" fill="#c5221f">AEX / Affinity Skid</text>

          <!-- Column Body -->
          <rect x="80" y="70" width="80" height="170" rx="6" fill="#ffffff" stroke="#c5221f" stroke-width="3"/>
          <rect x="84" y="100" width="72" height="120" fill="#fce8e6"/>

          <circle cx="100" cy="120" r="4" fill="#c5221f" opacity="0.5"/>
          <circle cx="130" cy="140" r="5" fill="#c5221f" opacity="0.5"/>
          <circle cx="140" cy="175" r="4" fill="#c5221f" opacity="0.5"/>
          <circle cx="110" cy="195" r="5" fill="#c5221f" opacity="0.5"/>

          <rect x="72" y="58" width="96" height="12" fill="#5f6368" rx="3"/>
          <rect x="72" y="240" width="96" height="12" fill="#5f6368" rx="3"/>

          <!-- Readable Text Box -->
          <rect x="90" y="125" width="60" height="42" rx="4" fill="#ffffff" opacity="0.9"/>
          <text x="120" y="142" font-family="Arial" font-size="10" font-weight="bold" text-anchor="middle" fill="#202124">Вихід: {downstream_recovery*100:.0f}%</text>
          <text x="120" y="158" font-family="Arial" font-size="11" font-weight="bold" text-anchor="middle" fill="#137333">% Full: {full_pct:.1f}%</text>
        </svg>
        """
    st.components.v1.html(svg_chrom, height=300)

  # 5. FILL/FINISH & QC (Ізолятор, Розлив, -80°C)
  with eq5:
    st.markdown("##### 5. Розлив & QC\n**(Isolator & -80°C Freezer)**")
    svg_fill = f"""
        <svg width="240" height="290" viewBox="0 0 240 290" xmlns="http://www.w3.org/2000/svg">
          <rect x="5" y="5" width="230" height="280" rx="12" fill="#f8f9fa" stroke="#8e24aa" stroke-width="2.5"/>
          <text x="120" y="32" font-family="Arial" font-size="13" font-weight="bold" text-anchor="middle" fill="#8e24aa">Class A Isolator</text>

          <!-- Isolator Box -->
          <rect x="35" y="55" width="170" height="105" rx="8" fill="#f3e8fd" stroke="#8e24aa" stroke-width="2"/>
          <text x="120" y="75" font-family="Arial" font-size="11" font-weight="bold" text-anchor="middle" fill="#202124">ddPCR / qPCR</text>
          <text x="120" y="93" font-family="Arial" font-size="11" font-weight="bold" text-anchor="middle" fill="#8e24aa">Score: {quality_score:.0f}/100</text>

          <!-- Glass Vials -->
          <rect x="55" y="110" width="18" height="32" rx="3" fill="#ffffff" stroke="#8e24aa" stroke-width="2"/>
          <rect x="88" y="110" width="18" height="32" rx="3" fill="#ffffff" stroke="#8e24aa" stroke-width="2"/>
          <rect x="121" y="110" width="18" height="32" rx="3" fill="#ffffff" stroke="#8e24aa" stroke-width="2"/>
          <rect x="154" y="110" width="18" height="32" rx="3" fill="#ffffff" stroke="#8e24aa" stroke-width="2"/>

          <rect x="58" y="124" width="12" height="16" fill="#8e24aa" opacity="0.6"/>
          <rect x="91" y="124" width="12" height="16" fill="#8e24aa" opacity="0.6"/>
          <rect x="124" y="124" width="12" height="16" fill="#8e24aa" opacity="0.6"/>
          <rect x="157" y="124" width="12" height="16" fill="#8e24aa" opacity="0.6"/>

          <!-- Freezer -->
          <rect x="45" y="180" width="150" height="70" rx="8" fill="#4a148c"/>
          <text x="120" y="220" font-family="Arial" font-size="14" font-weight="bold" text-anchor="middle" fill="#ffffff">Сховище -80°C</text>
        </svg>
        """
    st.components.v1.html(svg_fill, height=300)

  st.markdown("---")
  st.subheader("⚖️ Порівняльна характеристика масштабів виробництва")

  m1, m2, m3 = st.columns(3)

  with m1:
    st.markdown(
        """
        <div style="border:1px solid #dadce0; border-radius:8px; padding:15px; background:#ffffff;">
            <h3 style="text-align:center; color:#5f6368; margin-top:0;">1. LAB R&D SETUP</h3>
            <p style="text-align:center; color:#70757a; font-size:13px;">Масштаб 2 L – 10 L</p>
            <hr>
            <ul>
                <li><b>Обладнання:</b> Шейкерні колби, настільні реактори</li>
                <li><b>Гнучкість:</b> Дуже висока</li>
                <li><b>Основна мета:</b> Пошук параметрів трансфекції та тестування плазмід</li>
                <li><b>Вихід:</b> 10¹² – 10¹³ vg/лот</li>
            </ul>
        </div>
        """,
        unsafe_allow_html=True,
    )

  with m2:
    st.markdown(
        f"""
        <div style="border:2px solid #1a73e8; border-radius:8px; padding:15px; background:#f8f9fa;">
            <h3 style="text-align:center; color:#1a73e8; margin-top:0;">2. PILOT PRODUCTION</h3>
            <p style="text-align:center; color:#1a73e8; font-size:13px;">Масштаб 50 L – 200 L (Поточний: {scale_l} L)</p>
            <hr>
            <ul>
                <li><b>Обладнання:</b> Single-Use Bioreactor (SUB)</li>
                <li><b>Гнучкість:</b> Поміркована (дослідницькі лоти)</li>
                <li><b>Основна мета:</b> Доклінічні випробування (Toxicology studies)</li>
                <li><b>Вихід:</b> {purified_total_vg:.2e} vg/лот</li>
            </ul>
        </div>
        """,
        unsafe_allow_html=True,
    )

  with m3:
    st.markdown(
        """
        <div style="border:1px solid #dadce0; border-radius:8px; padding:15px; background:#ffffff;">
            <h3 style="text-align:center; color:#137333; margin-top:0;">3. COMMERCIAL GMP</h3>
            <p style="text-align:center; color:#70757a; font-size:13px;">Масштаб 1000 L – 2000 L</p>
            <hr>
            <ul>
                <li><b>Обладнання:</b> Промислові SUB-ферментери + Автоматизовані хроматографи</li>
                <li><b>Гнучкість:</b> Низька (Суворі GMP стандарти)</li>
                <li><b>Основна мета:</b> Клінічні фази та комерційний випуск</li>
                <li><b>Вихід:</b> >10¹⁶ vg/лот</li>
            </ul>
        </div>
        """,
        unsafe_allow_html=True,
    )

# --- ТАБ 3: КАРТКИ-ШПОРИ ---
with tab3:
  st.subheader(
      "📚 Протокольні картки-шпори біоінженерії (HEK293 Triple Transfection)"
  )

  col_card1, col_card2 = st.columns(2)

  with col_card1:
    st.info("""
        ### 🎯 1. Оптимізація трансфекції та комплексу
        
        * **Щільність суспензії HEK293:** Вводьте комплекс при **1.5–2.5 × 10⁶ cells/mL** із життєздатністю >95%.
        * **Дозування ДНК:** Концентрація загальної ДНК **1.0–1.5 мкг / 10⁶ клітин** дає оптимуми.
        * **Комплексоутворення DNA:PEI:** Співвідношення мас має бути **1:2 – 1:4 (w/w)**. Комплекс повинен формуватися строго **15–20 хвилин** перед додаванням. Понад 25 хвилин призводить до агломерації.
        * **Баланс плазмід:** Дотримуйтесь еквімолярного співвідношення `Transgene : Rep/Cap : Helper` **1:1:1**.
        """)

    st.success("""
        ### ❄️ 2. Температурний режим та Збір (Harvest)
        
        * **Гіпотермічний шифт (33°C):** Зниження температури з 37°C до 32–34°C після трансфекції покращує згортання білків, підвищує частку повних капсидів та зменшує метаболічний стрес.
        * **Час збору:** Оптимальне вікно — **48–72 hpt**. Збирайте біомасу до того, як життєздатність впаде нижче 70%.
        * **Серотипний розподіл:** AAV2 та AAV5 залишаються в клітинах (потрібен лізис), а AAV8 та AAV9 значною мірою вивільняються у середовище.
        """)

  with col_card2:
    st.warning("""
        ### 🧪 3. Лізис та Очищення (Downstream)
        
        * **Протоколи лізису:** Цикли заморожування-розморожування підходять для Small-Scale. Для біореакторів використовують м'який детергентний лізис або механічне руйнування.
        * **Очищення (Downstream):** 
          * *Без очищення* — високий вихід за геномами, але непридатний через HCDNA/HCP.
          * *Ультрацентрифугування (CsCl/Iodixanol)* — вихід ~40%, чистота ~85%.
          * *Афінна + Іонообмінна хроматографія (Affinity + IEX)* — для клінічного GMP-виробництва, чистота >99%, вихід ~65%.
        """)

    st.error("""
        ### 🏭 4. Чек-лист успішного Scale-Up (2000 L)
        
        * **P/V & kLa:** Масштабуйте процес за рівною питомою потужністю перемішування ($P/V$) та коефіцієнтом масопереносу кисню ($k_L a$).
        * **Inline-додавання комплексу:** У 2000 L застосовуйте автоматичну швидку подачу комплексу, щоб не вийти за межі 20 хвилин комплексоутворення.
        """)

# --- ТАБ 4: МАСШТАБУВАННЯ ТА ДИНАМІЧНЕ РЕЗЮМЕ ---
with tab4:
  st.subheader("⚖️ Порівняльний аналіз масштабування: 200 L ➔ 2000 L")

  case_data = {
      "Показник": [
          "Живі клітини при трансфекції",
          "Час збору",
          "Титр Crude Harvest (vg/mL)",
          "Частка повних капсидів (% full)",
          "Downstream-вихід (Purified)",
          "Критичне лімітування",
      ],
      "Пілотний біореактор (200 L)": [
          "2.0 x 10⁶ cells/mL",
          f"{harvest_hpt} hpt",
          f"{titer:.2e}",
          f"{full_pct:.1f}%",
          f"{purified_total_vg:.2e} vg",
          "Швидке однорідне змішування",
      ],
      "Промисловий біореактор (2000 L)": [
          "2.0 x 10⁶ cells/mL",
          f"{harvest_hpt} hpt",
          f"{titer * (0.78 / scale_penalty):.2e}",
          f"{max(5.0, full_pct - 3.0):.1f}%",
          f"{purified_total_vg * (0.78 / scale_penalty):.2e} vg",
          "Тривалість додавання комплексу, Shear stress, kLa",
      ],
  }
  st.dataframe(pd.DataFrame(case_data), use_container_width=True)

  st.markdown("---")
  st.subheader("💡 Динамічне інженерне резюме та рекомендації")

  if recommendations:
    st.warning(
        "**Система розрахувала наступні оптимізаційні кроки для покращення"
        " вашого лоту:**\n\n"
        + "\n\n".join(recommendations)
    )
  else:
    st.success(
        "🎉 **Ідеальні параметри біопроцесу!** Поточні налаштування забезпечують"
        " максимальну якість препарату (Quality Score), високу частку повних"
        " капсидів та високий чистий вихід вектора rAAV за стандартами GMP."
    )
