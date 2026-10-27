import pandas as pd
import streamlit as st

st.set_page_config(
    page_title="Навчальний симулятор rAAV", page_icon="🧬", layout="wide"
)

st.title("🧬 Навчальний симулятор напрацювання rAAV у біореакторі")
st.caption(
    "Платформа: HEK293 + транзієнтна потрійна плазмідна трансфекція. "
    "Модель показує компроміси між титром, якістю та параметрами масштабування."
)

# --- SIDEBAR: ВХІДНІ ПАРАМЕТРИ ---
st.sidebar.header("Параметри процесу")

cell_density = st.sidebar.slider(
    "Живі клітини при трансфекції, 10⁶ cells/mL",
    min_value=0.5,
    max_value=8.0,
    value=3.0,
    step=0.1,
)

dna_per_cell = st.sidebar.slider(
    "Загальна ДНК на клітину, ум. од.",
    min_value=0.2,
    max_value=3.0,
    value=1.0,
    step=0.1,
)

st.sidebar.subheader("Співвідношення плазмід (сума = 1.0)")
ratio_transgene = st.sidebar.slider(
    "Касета трансгена (GoI)", min_value=0.1, max_value=0.7, value=0.33, step=0.01
)
ratio_repcap = st.sidebar.slider(
    "Rep/Cap плазміда", min_value=0.1, max_value=0.7, value=0.33, step=0.01
)
ratio_helper = st.sidebar.slider(
    "Helper плазміда", min_value=0.1, max_value=0.7, value=0.34, step=0.01
)

harvest_hpt = st.sidebar.slider(
    "Час збору після трансфекції, hpt (годин)",
    min_value=24,
    max_value=120,
    value=60,
    step=1,
)

serotype = st.sidebar.selectbox(
    "Серотип rAAV", options=["AAV2", "AAV5", "AAV8", "AAV9"]
)

lysis_mode = st.sidebar.selectbox(
    "Стратегія лізису",
    options=[
        "Без лізису (збір супернатанту)",
        "Частковий лізис",
        "Повний лізис",
    ],
)

scale_l = st.sidebar.select_slider(
    "Масштаб біореактора, L", options=[2, 50, 200, 500, 1000, 2000], value=200
)

# --- РОЗРАХУНКИ ---

# Base titer (vg/mL)
base_titer = 3.0e10

# Ефект щільності клітин (оптимум ~3.0 x 10^6 cells/mL)
cell_effect = max(0.1, 1.0 - abs(cell_density - 3.0) / 8.0)

# Ефект кількості ДНК на клітину
dna_effect = max(
    0.1, 1.0 + 0.25 * (dna_per_cell - 1.0) - 0.20 * max(0.0, dna_per_cell - 1.5)
)

# Дисбаланс плазмід
plasmid_sum = ratio_transgene + ratio_repcap + ratio_helper
imbalance = abs(plasmid_sum - 1.0) + abs(ratio_transgene - 0.33) * 0.5
balance_effect = max(0.1, 1.0 - 1.5 * imbalance)

# Час збору
if harvest_hpt < 48:
  time_effect = 0.6 + (harvest_hpt - 24) / 24 * 0.4
elif harvest_hpt <= 72:
  time_effect = 1.0
else:
  time_effect = 1.0 - (harvest_hpt - 72) / 48 * 0.25

# Внутрішньоклітинна фракція за серотипами
intracellular_fraction = {
    "AAV2": 0.85,
    "AAV5": 0.60,
    "AAV8": 0.75,
    "AAV9": 0.70,
}[serotype]

lysis_release = {
    "Без лізису (збір супернатанту)": 0.0,
    "Частковий лізис": 0.55,
    "Повний лізис": 0.95,
}[lysis_mode]

available_fraction = (1.0 - intracellular_fraction) + (
    intracellular_fraction * lysis_release
)

# Штраф за масштаб (нелінійні втрати через mixing time, kLa, shear stress)
scale_penalty = {2: 1.00, 50: 0.98, 200: 0.95, 500: 0.90, 1000: 0.85, 2000: 0.78}[
    scale_l
]

# Кінцевий титр crude harvest (vg/mL)
titer = (
    base_titer
    * cell_effect
    * dna_effect
    * balance_effect
    * time_effect
    * available_fraction
    * scale_penalty
)

# Частка повних капсидів (%)
base_full = 0.28
full_fraction = (
    base_full
    - 0.08 * max(0.0, dna_per_cell - 1.0)
    - 0.15 * imbalance
    - 0.05 * max(0.0, harvest_hpt - 72) / 24.0
    - 0.03 * max(0.0, cell_density - 5.0)
)
full_pct = max(5.0, min(45.0, full_fraction * 100))

# Домішки
impurity_load = (
    1.0
    + 2.0 * lysis_release
    + 0.5 * max(0.0, harvest_hpt - 72) / 24.0
    + 0.3 * max(0.0, cell_density - 4.0)
)

# Індекс якості (0-100)
quality_score = (
    0.50 * (full_pct / 30.0 * 100)
    + 0.30 * max(0.0, (1.0 - (impurity_load - 1.0) / 3.0) * 100)
    + 0.20 * min(100.0, (titer / 3.0e10) * 100)
)
quality_score = max(0.0, min(100.0, quality_score))

# Загальний вихід (vg)
harvest_volume_ml = scale_l * 1000.0 * 0.9
total_vg = titer * harvest_volume_ml

# --- РЕЗУЛЬТАТИ ---
st.subheader("📊 Прогнозовані показники процесу")
c1, c2, c3, c4 = st.columns(4)
c1.metric("Геномний титр (Crude)", f"{titer:.2e} vg/mL")
c2.metric("Повний вихід", f"{total_vg:.2e} vg")
c3.metric("Частка повних капсидів", f"{full_pct:.1f} %")
c4.metric("Індекс якості лоту", f"{quality_score:.0f} / 100")

# --- ЗАСТЕРЕЖЕННЯ ТА КОМПРОМІСИ ---
st.subheader("⚠️ Аналіз компромісів та застереження")

if abs(plasmid_sum - 1.0) > 0.02:
  st.error(
      f"Сума часток плазмід має дорівнювати 1.00. Наразі сума:"
      f" {plasmid_sum:.2f}."
  )

if titer > 3.0e10 and full_pct < 20.0:
  st.warning(
      "⚡ **Титр високий, але частка повних капсидів знизилася!** "
      "Пам'ятайте, що геномний титр (vg/mL) не дорівнює інфекційному титру."
  )

if dna_per_cell > 1.8:
  st.warning(
      "🧬 **Надлишок ДНК на клітину:** підвищує кількість порожніх або"
      " часткових капсидів і створює навантаження на очищення."
  )

if harvest_hpt > 84:
  st.warning(
      "⏱️ **Пізній збір (>84 hpt):** пришвидшує деградацію капсиду та збільшує"
      " кількість домішок хазяйських білків і ДНК (HCP/HCDNA)."
  )

if lysis_mode == "Повний лізис":
  st.warning(
      "🧪 **Повний лізис:** максимізує вивільнення вектора, але різко"
      " збільшує вміст хазяйської ДНК і білка, які потрібно видаляти на"
      " стадії downstream."
  )

if scale_l >= 1000:
  st.info(
      "🏭 **Масштаб 1000–2000 L:** час додавання трансфекційного комплексу,"
      " неоднорідність змішування та кисневого масопереносу (kLa) знижують"
      " питому продуктивність."
  )

# --- КЕЙС-СТАДІ МАСШТАБУВАННЯ ---
st.subheader("⚖️ Кейс-стаді масштабування: 200 L ➔ 2000 L")
case_data = {
    "Показник": [
        "Живі клітини при трансфекції",
        "Час збору",
        "Прогнозований титр (vg/mL)",
        "Частка повних капсидів (% full)",
        "Навантаження домішками",
        "Обмежувальний фактор",
    ],
    "Лабораторний / Пілот (200 L)": [
        "3.0 x 10⁶ cells/mL",
        f"{harvest_hpt} hpt",
        f"{titer:.2e}",
        f"{full_pct:.1f}%",
        "Базове",
        "Швидке змішування, стабільний DO",
    ],
    "Промисловий (2000 L)": [
        "3.0 x 10⁶ cells/mL",
        f"{harvest_hpt} hpt",
        f"{titer * (0.78 / scale_penalty):.2e}",
        f"{max(5.0, full_pct - 3.0):.1f}%",
        "Підвищене",
        "Тривалість додавання комплексу, зсувне напруження, kLa",
    ],
}
st.dataframe(pd.DataFrame(case_data), use_container_width=True)

st.success(
    "💡 **Головний висновок scale-up:** критерієм успіху є збереження якості та"
    " продуктивності на літр, а не досягнення максимального vg/mL за будь-яку"
    " ціну."
)
